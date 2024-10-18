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
    def __init__(self, feature_label, entity_type, id_type, file_path, selected_column, parent=None):
        super().__init__(parent)

        self.feature_label = feature_label
        self.entity_type = entity_type
        self.id_type = id_type
        self.file_path = file_path
        self.selected_column = selected_column

        # Layout for the row
        self.layout = QVBoxLayout(self)

        # Display feature label, entity type and file path
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
            # Start the process in a separate thread
            self.parentWidget().start_processing(self, process_func)
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
        self.layout = QVBoxLayout(self)
        self.process_rows = []
        self.selected_entities = {}

        # Connect the signal to the slot that shows the dialog
        self.show_dialog_signal.connect(self.show_multi_dialog)

        # Add loading label
        self.loading_label = QLabel("Processing data...")
        self.layout.addWidget(self.loading_label)
        self.animation_chars = cycle(["⢿", "⣻", "⣽", "⣾", "⣷", "⣯", "⣟", "⡿"])
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_loading_animation)
        self.loading_label.hide()  # Initially hide the label

        # Create a ProcessRow for each entry in read_info_list if it's provided
        if read_info_list:
            for feature_label, entity_type, id_type, file_path, selected_column in read_info_list:
                row = ProcessRow(feature_label, entity_type, id_type, file_path, selected_column)
                self.process_rows.append(row)
                self.layout.addWidget(row)

        # Add cache path section
        self.add_cache_path_section()

        # Add some spacing at the bottom
        self.layout.addStretch()

        # Create thread pool
        self.threadpool = QThreadPool()

    def update_loading_animation(self):
        """Update the loading label with the next animation character."""
        char = next(self.animation_chars)
        self.loading_label.setText(f"Processing data... {char}")

    def start_processing(self, row, process_func):
        """Start processing data in the background."""
        self.loading_label.show()  # Show the loading label
        self.timer.start(100)  # Start the animation timer
        
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
            # none matcher func
            worker = Worker(
                process_func,
                row.entity_type,
                row.id_type,
                row.file_path,
                row.selected_column,
                row.feature_label
            )

        worker.signals.finished.connect(self.on_processing_complete)
        worker.signals.error.connect(self.on_processing_error)
        self.threadpool.start(worker)

    def emit_show_dialog_signal(self, all_topk_phenotypes, entity_type, id_type, file_path, selected_column, feature_label, phenotype):
        """Emit signal to show dialog in the main thread."""
        self.entity_type = entity_type
        self.id_type = id_type
        self.file_path = file_path
        self.selected_column = selected_column
        self.feature_label = feature_label
        self.phenotype = phenotype

        self.show_dialog_signal.emit(all_topk_phenotypes)

    def show_multi_dialog(self, all_topk_entities):
        """Display the multi-selection dialog in the main thread."""
        print("Preparing to show multi-selection dialog for all entities.")

        # Pop up a dialog for multi-selection
        dialog = MultiEntitySelectionDialog(all_topk_entities)
        if dialog.exec_() == QDialog.Accepted:
            selected_entities = dialog.get_selected_entities()
            if selected_entities:
                self.selected_entities = selected_entities
                print(f"User selections: {self.selected_entities}")
                # Continue processing after selection
                self.continue_processing_after_selection()
        else:
            print("No selections made")

    def continue_processing_after_selection(self):
        """Pass the selected entities back for further processing."""
        process_entities_file(
            self.entity_type,
            self.id_type,
            self.feature_label,
            self.phenotype,
            self.selected_column,
            self.selected_entities
        )

        self.timer.stop()
        self.loading_label.setText("Processing complete.")
        QTimer.singleShot(2000, self.loading_label.hide)

    def on_processing_complete(self):
        """Handle processing completion."""
        self.timer.stop()  # Stop the animation timer
        self.loading_label.setText("Processing complete.")
        QTimer.singleShot(2000, self.loading_label.hide)  # Hide after 2 seconds

    def on_processing_error(self, error):
        """Handle processing error."""
        self.timer.stop()  # Stop the animation timer
        exctype, value, traceback_str = error
        QMessageBox.critical(self, "Error", f"An error occurred: {value}\n{traceback_str}")
        self.loading_label.hide()

    def add_cache_path_section(self):
        """Add a section at the bottom for Cache path and controls."""
        cache_layout = QVBoxLayout()

        # Cache path label
        cache_label = QLabel("Cache path:")
        cache_layout.addWidget(cache_label)

        # Create the horizontal layout for the input and buttons
        path_layout = QHBoxLayout()

        # Path input field
        self.cache_path_input = QLineEdit("./cache")  # Default to ./cache
        path_layout.addWidget(self.cache_path_input)

        # Browse button
        self.browse_button = QPushButton("Browse")
        self.browse_button.clicked.connect(self.browse_cache_path)
        path_layout.addWidget(self.browse_button)

        # Process all button
        self.process_all_button = QPushButton("Process All")
        self.process_all_button.clicked.connect(self.process_all_data)
        path_layout.addWidget(self.process_all_button)

        # Add the path layout to the cache layout
        cache_layout.addLayout(path_layout)

        # Add some spacing at the bottom
        cache_layout.addStretch()

        # Add the cache layout to the main layout
        self.layout.addLayout(cache_layout)

    def browse_cache_path(self):
        """Open a file dialog to select the cache directory."""
        directory = QFileDialog.getExistingDirectory(self, "Select Cache Directory", "./")
        if directory:
            self.cache_path_input.setText(directory)

    def process_all_data(self):
        """Handle processing all data from the cache path in a separate thread."""
        cache_path = self.cache_path_input.text()
        
        # Start the process in a separate thread
        worker = Worker(process_all_data, cache_path)
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
        self.table_widget.setColumnWidth(1, 400)

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