# process_tab.py

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from itertools import cycle
from data_process import *

from entity_process.gene_process import *
from entity_process.transcript_process import *
from entity_process.protein_process import *
from entity_process.promoter_process import *
from entity_process.drug_process import *
from entity_process.disease_process import *
from entity_process.phenotype_process import *
from entity_process.soft_match_process import *

import sys
import traceback

class WorkerSignals(QObject):
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)

class Worker(QRunnable):
    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
        except Exception as e:
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()

class ProcessRow(QWidget):
    def __init__(self, feature_label, entity_type, id_type, file_path, selected_column, process_tab_ref, parent=None):
        super().__init__(parent)

        self.feature_label = feature_label
        self.entity_type = entity_type
        self.id_type = id_type
        self.file_path = file_path
        self.selected_column = selected_column
        self.process_tab_ref = process_tab_ref  # Reference to ProcessTab

        # Layout for the row
        self.layout = QVBoxLayout(self)

        # Display feature label, entity type, and file path
        entity_layout = QHBoxLayout()
        feature_label_display = QLabel(f"Feature Label: {feature_label}")
        entity_label = QLabel(f"Entity Type: {entity_type}")
        self.process_button = QPushButton("Process")
        self.process_button.setFixedWidth(200)

        # Add widgets to the layout
        entity_layout.addWidget(feature_label_display)
        entity_layout.addWidget(entity_label)
        entity_layout.addWidget(self.process_button)

        # File path display
        path_label = QLabel(file_path)
        path_label.setStyleSheet("color: gray; font-size: 12px;")

        self.layout.addLayout(entity_layout)
        self.layout.addWidget(path_label)

        # Connect process button to the appropriate function
        self.process_button.clicked.connect(self.process_data)

    def process_data(self):
        """Call the appropriate process function based on the entity type."""
        process_functions = {
            "Gene": process_gene,
            "Transcript": process_transcript,
            "Protein": process_protein,
            "Promoter": process_promoter,
            "Drug": process_drug,
            "Disease": process_disease,
            "Phenotype": process_phenotype,
        }

        process_func = process_functions.get(self.entity_type)
        if process_func:
            # Start the process using the reference to ProcessTab
            self.process_tab_ref.start_processing(self, process_func)
        else:
            QMessageBox.critical(self, "Error", f"No processing function found for entity type: {self.entity_type}")
    
class ProcessTab(QWidget):

    show_dialog_signal = pyqtSignal(dict)

    def __init__(self, read_info_list=None, matcher=None, phenotype_embeddings=None, disease_embeddings=None, drug_embeddings=None, parent=None):
        super().__init__(parent)
        self.matcher = matcher
        self.phenotype_embeddings = phenotype_embeddings
        self.disease_embeddings = disease_embeddings
        self.drug_embeddings = drug_embeddings

        # Create a main layout for the tab
        main_layout = QVBoxLayout(self)

        # Create a QSplitter to divide the interface into two sections
        splitter = QSplitter(Qt.Vertical)  # Vertical splitter for top and bottom

        # Upper section for processing rows (add scroll area for multiple rows)
        upper_widget = QWidget()
        upper_layout = QVBoxLayout(upper_widget)
        self.process_rows = []

        if read_info_list:
            for feature_label, entity_type, id_type, file_path, selected_column in read_info_list:
                # Pass self (reference to ProcessTab) to each ProcessRow
                row = ProcessRow(feature_label, entity_type, id_type, file_path, selected_column, process_tab_ref=self)
                self.process_rows.append(row)
                upper_layout.addWidget(row)

        # Add some spacing and stretch at the bottom of the upper section
        upper_layout.addStretch()

        # Scroll area for the upper section (in case there are many rows)
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(upper_widget)

        # Add scroll area to the upper part of the splitter
        splitter.addWidget(scroll_area)

        # Lower section for clinical data and cache path selection
        lower_widget = QWidget()
        lower_layout = QVBoxLayout(lower_widget)

        # Add clinical data section
        self.add_clinical_data_section(lower_layout)

        # Add cache path section
        self.add_cache_path_section(lower_layout)

        # Add some spacing at the bottom of the lower section
        lower_layout.addStretch()

        # Add lower widget to the splitter
        splitter.addWidget(lower_widget)

        # Adjust splitter size ratio
        splitter.setSizes([500, 150])  # Set initial sizes for top and bottom sections

        # Add the splitter to the main layout
        main_layout.addWidget(splitter)

        # Create thread pool
        self.threadpool = QThreadPool()

        # Add loading label to indicate when processing is in progress
        self.loading_label = QLabel("Processing data...")  # Define the loading label
        self.loading_label.setAlignment(Qt.AlignCenter)
        self.loading_label.hide()  # Initially hide the label
        main_layout.addWidget(self.loading_label)  # Add it to the main layout (central position)

        # Create animation characters for the loading animation
        self.animation_chars = cycle(["⢿", "⣻", "⣽", "⣾", "⣷", "⣯", "⣟", "⡿"])
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_loading_animation)

    def update_loading_animation(self):
        """Update the loading label with the next animation character."""
        char = next(self.animation_chars)
        self.loading_label.setText(f"Processing data... {char}")

    def start_processing(self, row, process_func):
        """Start processing data in the background (upper section)."""
        self.show_loading_animation()  # Show the loading label

        # Check the entity type and pass the appropriate embeddings
        if row.entity_type == "Phenotype":
            embeddings = self.phenotype_embeddings
        elif row.entity_type == "Disease":
            embeddings = self.disease_embeddings
        elif row.entity_type == "Drug":
            embeddings = self.drug_embeddings
        else:
            embeddings = None  # Handle other cases if necessary

        # Check if the id_type contains "name" to determine the process function
        if "name" in row.id_type.lower():
            # If id_type contains "name", use process_entities_file
            worker = Worker(
                process_entities,  # Call process_entities
                row.entity_type,
                row.id_type,
                row.file_path,
                row.selected_column,
                row.feature_label,
                matcher=self.matcher,
                embeddings=embeddings,
                selection_callback=self.emit_show_dialog_signal
            )
        else:
            # Use normal processing function
            worker = Worker(
                process_func,
                row.entity_type,
                row.id_type,
                row.file_path,
                row.selected_column,
                row.feature_label
            )

        worker.signals.finished.connect(self.on_processing_complete)
        worker.signals.error.connect(self.on_processing_error)  # Connect error handling
        self.threadpool.start(worker)

    def show_loading_animation(self):
        """Show and start the loading animation."""
        self.loading_label.show()
        self.timer.start(100)

    def on_processing_complete(self):
        """Handle processing completion."""
        self.timer.stop()
        self.loading_label.setText("Processing complete.")
        QTimer.singleShot(2000, self.loading_label.hide)  # Hide after 2 seconds

    def on_processing_error(self, error):
        """Handle processing error."""
        self.timer.stop()
        exctype, value, traceback_str = error
        QMessageBox.critical(self, "Error", f"An error occurred: {value}\n{traceback_str}")
        self.loading_label.hide()

    def add_clinical_data_section(self, layout):
        """Add a section at the bottom for Clinical data and controls."""
        clinical_layout = QVBoxLayout()

        # Clinical data label
        clinical_label = QLabel("Clinical data:")
        clinical_layout.addWidget(clinical_label)

        # Create the horizontal layout for the input and buttons
        clinical_data_path_layout = QHBoxLayout()

        # Data input field
        self.clinical_data_path_input = QLineEdit("./test_data/clinical_data.csv")
        
        # Browse button
        self.browse_button = QPushButton("Browse")
        self.browse_button.clicked.connect(self.browse_clinical_data_path)

        # Add input field and browse button to the horizontal layout
        clinical_data_path_layout.addWidget(self.clinical_data_path_input)
        clinical_data_path_layout.addWidget(self.browse_button)

        # Add the horizontal layout to the vertical layout
        clinical_layout.addLayout(clinical_data_path_layout)

        # Add the clinical layout to the passed layout
        layout.addLayout(clinical_layout)

    def add_cache_path_section(self, layout):
        """Add a section at the bottom for Cache path and controls."""
        cache_layout = QVBoxLayout()

        # Cache path label
        cache_label = QLabel("Cache path:")
        cache_layout.addWidget(cache_label)

        # Create the horizontal layout for the input and buttons
        cache_path_layout = QHBoxLayout()

        # Path input field
        self.cache_path_input = QLineEdit("./cache")  # Default to ./cache
        cache_path_layout.addWidget(self.cache_path_input)

        # Browse button
        self.browse_button = QPushButton("Browse")
        self.browse_button.clicked.connect(self.browse_cache_path)
        cache_path_layout.addWidget(self.browse_button)

        # Process all button
        self.process_all_button = QPushButton("Process All")
        self.process_all_button.clicked.connect(self.process_all_data)
        cache_path_layout.addWidget(self.process_all_button)

        # Add the path layout to the cache layout
        cache_layout.addLayout(cache_path_layout)

        # Add some spacing at the bottom
        cache_layout.addStretch()

        # Add the cache layout to the passed layout
        layout.addLayout(cache_layout)

    def browse_cache_path(self):
        """Open a file dialog to select the cache directory."""
        directory = QFileDialog.getExistingDirectory(self, "Select Cache Directory", "./")
        if directory:
            self.cache_path_input.setText(directory)

    def browse_clinical_data_path(self):
        """Open a file dialog to select a CSV file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Select Clinical Data File", 
            "./", 
            "CSV Files (*.csv)"  # Filter to show only .csv files
        )
        if file_path:
            self.clinical_data_path_input.setText(file_path)

    def process_all_data(self):
        """Handle processing all data from the cache path in a separate thread."""
        cache_path = self.cache_path_input.text()  # Get the cache path from input field
        clinical_data_path = self.clinical_data_path_input.text()  # Get the clinical data path from input field

        # Start the process in a separate thread and pass both cache_path and clinical_data_path
        worker = Worker(process_and_merge_data, cache_path, clinical_data_path)
        worker.signals.finished.connect(self.on_processing_complete)
        worker.signals.error.connect(self.on_processing_error)
        self.threadpool.start(worker)



class MultiEntitySelectionDialog(QDialog):
    def __init__(self, all_topk_entities, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Correct Entities")

        self.resize(1000, 800)

        self.selected_entities = {}

        # Layout for the dialog
        self.layout = QVBoxLayout(self)

        # Instruction label
        instruction_label = QLabel("Select the correct match for each entity:")
        self.layout.addWidget(instruction_label)

        # Table to display entities and options
        self.table_widget = QTableWidget(self)
        self.table_widget.setRowCount(len(all_topk_entities))
        self.table_widget.setColumnCount(2)  # Column 0: Entity, Column 1: Options
        self.table_widget.setHorizontalHeaderLabels(["Entity", "Select Matching Term"])

        # Set column widths
        self.table_widget.setColumnWidth(0, 300)
        self.table_widget.setColumnWidth(1, 500)

        # Populate the table with entity options
        for row_idx, (entity_value, topk_entities) in enumerate(all_topk_entities.items()):
            # Entity column
            entity_item = QTableWidgetItem(entity_value)
            self.table_widget.setItem(row_idx, 0, entity_item)

            # ComboBox for selecting matching term
            combo_box = QComboBox()
            for med_id, hpo_term in topk_entities:
                combo_box.addItem(f"{med_id}, {hpo_term}", (med_id, hpo_term))
            self.table_widget.setCellWidget(row_idx, 1, combo_box)

        self.layout.addWidget(self.table_widget)

        # OK and Cancel buttons
        button_layout = QHBoxLayout()
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")
        button_layout.addWidget(self.ok_button)
        button_layout.addWidget(self.cancel_button)
        self.layout.addLayout(button_layout)

        # Connect buttons to actions
        self.ok_button.clicked.connect(self.accept)
        self.cancel_button.clicked.connect(self.reject)

    def get_selected_entities(self):
        """Return the selected entities as a dictionary."""
        selected_entities = {}
        for row_idx in range(self.table_widget.rowCount()):
            entity_value = self.table_widget.item(row_idx, 0).text()
            combo_box = self.table_widget.cellWidget(row_idx, 1)
            selected_entity = combo_box.currentData()  # Get selected med_id, hpo_term tuple
            selected_entities[entity_value] = [selected_entity]
        return selected_entities