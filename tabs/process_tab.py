import os
import sys
import traceback
import shutil
import pandas as pd
from itertools import cycle
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from data_pipeline.data_process import *
from data_pipeline.entity_process.gene_process import *
from data_pipeline.entity_process.transcript_process import *
from data_pipeline.entity_process.protein_process import *
from data_pipeline.entity_process.promoter_process import *
from data_pipeline.entity_process.drug_process import *
from data_pipeline.entity_process.disease_process import *
from data_pipeline.entity_process.phenotype_process import *
from data_pipeline.entity_process.soft_match_process import *
from data_pipeline.path_manager import DatabasePathManager

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

        # Add process button to the parent tab's button list
        if self.process_tab_ref:
            self.process_tab_ref.process_buttons.append(self.process_button)

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
            # Start the process, passing the database_path
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
        self.process_rows = []
        self.selected_entities = {}
        self.file_order = []

        self.is_animation_running = False
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_loading_animation)

        # manage the buttons
        self.process_buttons = []

        # Load database path
        database_path = DatabasePathManager.get_database_path()
        if not database_path or not os.path.exists(database_path):
            QMessageBox.warning(self, "Warning", "Database path is not set or invalid. Please select a valid path.")
        else:
            print(f"Database path loaded: {database_path}")  # Debugging

        # Connect the signal to the slot that shows the dialog
        self.show_dialog_signal.connect(self.show_multi_dialog)

        # Create main layout for the tab
        main_layout = QVBoxLayout(self)

        # Initialize UI
        self.init_ui(main_layout, read_info_list)

        # Create thread pool
        self.threadpool = QThreadPool()

        # Add loading label
        self.loading_label = QLabel("Processing data...")
        self.loading_label.setAlignment(Qt.AlignLeft | Qt.AlignBottom)
        self.loading_label.setStyleSheet("font-weight: bold; color: black;")
        self.loading_label.hide()
        main_layout.addWidget(self.loading_label, alignment=Qt.AlignLeft)

        self.animation_chars = cycle(["⢿", "⣻", "⣽", "⣾", "⣷", "⣯", "⣟", "⡿"])
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_loading_animation)

    def validate_path(self):
        """Validate the database path."""
        path = DatabasePathManager.get_database_path()
        if not path or not os.path.exists(path):
            QMessageBox.critical(self, "Error", "Database path is not set or invalid. Please select a valid path.")
            return False
        return True

    def init_ui(self, main_layout, read_info_list):
        """Initialize the UI layout for the tab."""
        splitter = QSplitter(Qt.Vertical)

        # Upper section for processing rows
        upper_widget = QWidget()
        upper_layout = QVBoxLayout(upper_widget)

        upper_control_layout = QHBoxLayout()
        upper_control_layout.addStretch()
        upper_layout.addLayout(upper_control_layout)

        if read_info_list:
            for feature_label, entity_type, id_type, file_path, selected_column in read_info_list:
                row = ProcessRow(feature_label, entity_type, id_type, file_path, selected_column, process_tab_ref=self)
                self.process_rows.append(row)
                upper_layout.addWidget(row)

        upper_layout.addStretch()
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(upper_widget)
        splitter.addWidget(scroll_area)

        # Lower section for cache path
        lower_widget = QWidget()
        lower_layout = QVBoxLayout(lower_widget)

        order_layout = QHBoxLayout()
        self.entity_order_button = QPushButton("Entity Ordering")
        self.entity_order_button.clicked.connect(self.open_entity_order_dialog)
        self.entity_order_label = QLabel("No order selected")
        order_layout.addWidget(self.entity_order_button)
        order_layout.addWidget(self.entity_order_label)
        lower_layout.addLayout(order_layout)

        self.add_cache_path_section(lower_layout)
        lower_layout.addStretch()
        splitter.addWidget(lower_widget)
        splitter.setSizes([500, 150])

        main_layout.addWidget(splitter)

    def update_loading_animation(self):
        """Update the loading label with the next animation character."""
        char = next(self.animation_chars)
        current_text = self.loading_label.text().rstrip(" ⢿⣻⣽⣾⣷⣯⣟⡿")
        self.loading_label.setText(f"{current_text} {char}")

    def start_processing(self, row, process_func):
        """Start processing data."""
        if not self.validate_path():
            return  # Exit if path is invalid

        database_path = DatabasePathManager.get_database_path()
        task_name = f"Processing {row.entity_type} Data"

        # Show loading animation
        self.show_loading_animation(task_name)

        # Disable buttons during processing
        self.disable_process_buttons()

        embeddings = None
        if row.entity_type == "Phenotype":
            embeddings = self.phenotype_embeddings
        elif row.entity_type == "Disease":
            embeddings = self.disease_embeddings
        elif row.entity_type == "Drug":
            embeddings = self.drug_embeddings

        if "name" in row.id_type.lower():
            worker = Worker(
                process_entities,
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
            worker = Worker(
                process_func,
                row.entity_type,
                row.id_type,
                row.file_path,
                row.selected_column,
                row.feature_label,
                database_path=database_path
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
        dialog = MultiEntitySelectionDialog(all_topk_entities)
        if dialog.exec_() == QDialog.Accepted:
            selected_entities = dialog.get_selected_entities()
            if selected_entities:
                self.selected_entities = selected_entities
                self.continue_processing_after_selection()

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
        # QTimer.singleShot(2000, self.loading_label.hide)

    def show_loading_animation(self, task_name="Processing"):
        """Show the loading animation and disable buttons."""
        self.loading_label.setText(f"{task_name}...")
        self.loading_label.show()

        # Start the animation timer if not already running
        if not self.is_animation_running:
            self.is_animation_running = True
            self.timer.start(100)

    def on_processing_complete(self):
        """Handle processing completion."""
        self.timer.stop()
        self.is_animation_running = False
        self.loading_label.setText("Processing complete.")
        # QTimer.singleShot(2000, self.loading_label.hide)
        self.enable_process_buttons()

    def on_processing_error(self, error):
        """Handle processing errors."""
        self.timer.stop()
        self.is_animation_running = False
        self.loading_label.hide()

        exctype, value, traceback_str = error
        QMessageBox.critical(self, "Error", f"An error occurred: {value}\n{traceback_str}")

        # Re-enable buttons after processing
        self.enable_process_buttons()

    def open_entity_order_dialog(self):
        """Open the Entity Order Arrangement dialog and handle the result."""
        dialog = FileOrderDialog()  # Assume this gets files from ./cache
        if dialog.exec_() == QDialog.Accepted:
            self.file_order = dialog.ordered_files  # Store the file order to self.file_order
            # Display selected order in the label with arrows
            if self.file_order:
                self.entity_order_label.setText(" → ".join(self.file_order))
            else:
                self.entity_order_label.setText("No order selected")
        else:
            self.entity_order_label.setText("No order selected")

    def get_file_order(self):
        """Return the file order for other functions to use."""
        return self.file_order

    def add_cache_path_section(self, layout):
        """Add a section at the bottom for Cache path and controls."""
        cache_layout = QVBoxLayout()

        # Cache path label
        cache_label = QLabel("Cache path (default: ./cache)")
        cache_layout.addWidget(cache_label)

        # Create the horizontal layout for the input and browse button
        cache_path_layout = QHBoxLayout()

        # Path input field
        self.cache_path_input = QLineEdit("./cache")  # Default to ./cache
        cache_path_layout.addWidget(self.cache_path_input)

        # Browse button (excluded from disable logic)
        self.browse_button = QPushButton("Browse")
        self.browse_button.clicked.connect(self.browse_cache_path)
        cache_path_layout.addWidget(self.browse_button)

        # Add the path layout to the cache layout
        cache_layout.addLayout(cache_path_layout)

        # Finalize button
        self.finalize_button = QPushButton("Finalize")
        self.finalize_button.clicked.connect(self.finalize_data)
        self.finalize_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)  # Expand horizontally
        self.process_buttons.append(self.finalize_button)  # Save button reference
        cache_layout.addWidget(self.finalize_button)

        # Add some spacing at the bottom
        cache_layout.addStretch()

        # Add the cache layout to the passed layout
        layout.addLayout(cache_layout)

    def browse_cache_path(self):
        """Open a file dialog to select the cache directory."""
        directory = QFileDialog.getExistingDirectory(self, "Select Cache Directory", "./")
        if directory:
            self.cache_path_input.setText(directory)

    def finalize_data(self):
        """Handle processing all data in two steps."""
        cache_path = self.cache_path_input.text()
        file_order = self.get_file_order()

        # Check if file_order is empty and prompt the user
        if not file_order:
            QMessageBox.warning(self, "Entity Order Required", "Please arrange the Entity Order before proceeding.")
            return  # Exit the function if file_order is empty

        # Disable buttons during processing
        self.disable_process_buttons()

        self.show_loading_animation("Finalizing data")

        # Merge and prepare mapping
        worker = Worker(self.merge_and_prepare_mapping, cache_path, file_order)
        worker.signals.result.connect(self.on_mapping_prepared)
        worker.signals.error.connect(self.on_processing_error)
        worker.signals.finished.connect(self.enable_process_buttons)  # Re-enable buttons after processing
        self.threadpool.start(worker)

    def merge_and_prepare_mapping(self, cache_path, file_order):
        """First step: merge data and prepare entity index mapping."""
        entity_index_id_mapping_df, merged_data, processed_data_path = merge_data_and_generate_entity_mapping(
            cache_path, file_order
        )
        return entity_index_id_mapping_df, merged_data, processed_data_path

    def on_mapping_prepared(self, result):
        """Callback after merging and entity mapping generation."""
        database_path = DatabasePathManager.get_database_path()
        entity_index_id_mapping_df, merged_data, processed_data_path = result

        if entity_index_id_mapping_df is None:
            QMessageBox.critical(self, "Error", "Failed to prepare entity mapping.")
            self.on_processing_complete()
            return

        # Filter and save edge data
        worker = Worker(self.filter_edge_data_types, database_path, entity_index_id_mapping_df)
        worker.signals.result.connect(self.on_unique_types_fetched)
        worker.signals.error.connect(self.on_processing_error)
        self.threadpool.start(worker)

    def filter_edge_data_types(self, database_path, entity_index_id_mapping_df):
        """Filter edge data and return unique Type values."""
        filtered_edge_data, unique_types = filter_and_save_edge_data(database_path, entity_index_id_mapping_df)
        return filtered_edge_data, unique_types

    def on_unique_types_fetched(self, result):
        """Show dialog to select edge data types based on unique values."""
        filtered_edge_data, unique_types = result

        # Show dialog to select edge types
        dialog = TypeSelectionDialog(unique_types)
        if dialog.exec_() == QDialog.Accepted:
            selected_types = dialog.get_selected_types()
            
            # Filter edge data based on selected types
            self.process_selected_edge_types(filtered_edge_data, selected_types)
        else:
            QMessageBox.warning(self, "Cancelled", "Edge data processing was cancelled.")
            self.on_processing_complete()

    def process_selected_edge_types(self, filtered_edge_data, selected_types):
        """Filter edge data by selected types and save results."""
        worker = Worker(self.final_edge_processing, filtered_edge_data, selected_types)
        worker.signals.finished.connect(self.on_processing_complete)
        worker.signals.error.connect(self.on_processing_error)
        self.threadpool.start(worker)

    def final_edge_processing(self, filtered_edge_data, selected_types):
        """Final processing of edge data based on selected types."""
        # Use the selected types to process the edge data
        entity_index_id_mapping_df = pd.read_csv('./cache/processed_data/entity_index_id_mapping.csv')
        process_edge_data_with_selected_types(filtered_edge_data, selected_types, entity_index_id_mapping_df, './cache/processed_data')

    def disable_process_buttons(self):
        """Disable all process-related buttons."""
        for button in self.process_buttons:
            button.setEnabled(False)

    def enable_process_buttons(self):
        """Enable all process-related buttons."""
        for button in self.process_buttons:
            button.setEnabled(True)

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
        self.setWindowIcon(QIcon("assets/icons/logo.png"))

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
        self.cancel_button = QPushButton("Cancel") # left button
        self.ok_button = QPushButton("OK") # right button
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
    
class FileOrderDialog(QDialog):
    def __init__(self, cache_dir="./cache", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Arrange Entity Order")
        self.resize(400, 300)

        # Main layout
        layout = QVBoxLayout(self)

        # Create a QListWidget to display files
        self.file_list_widget = QListWidget()
        self.file_list_widget.setDragDropMode(QListWidget.InternalMove)  # Enable drag and drop to reorder
        layout.addWidget(self.file_list_widget)
        self.setWindowIcon(QIcon("assets/icons/logo.png"))

        # Populate the QListWidget with files from cache directory
        self.populate_file_list(cache_dir)

        # Add a button to confirm the order
        save_button = QPushButton("Confirm Order")
        save_button.clicked.connect(self.confirm_order)
        layout.addWidget(save_button)

        self.ordered_files = []  # To store the final order of files

    def populate_file_list(self, cache_dir):
        """Populate QListWidget with files from the cache directory."""
        if not os.path.exists(cache_dir):
            QMessageBox.warning(self, "Error", f"The directory '{cache_dir}' does not exist.")
            return

        files = os.listdir(cache_dir)
        files = [f for f in files if os.path.isfile(os.path.join(cache_dir, f))]  # Only add files, not directories

        if not files:
            QMessageBox.warning(self, "Error", f"No files found in the directory '{cache_dir}'.")
            return

        for file_name in files:
            item = QListWidgetItem(file_name)
            self.file_list_widget.addItem(item)

    def confirm_order(self):
        """Save the file order after the user has reordered the list."""
        self.ordered_files = []
        for i in range(self.file_list_widget.count()):
            item = self.file_list_widget.item(i)
            file_name_without_ext = os.path.splitext(item.text())[0]  # Get file name without extension
            self.ordered_files.append(file_name_without_ext)
        
        self.accept()  # Close the dialog and return the result


class TypeSelectionDialog(QDialog):
    def __init__(self, unique_types, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Edge Types")
        self.setWindowIcon(QIcon("assets/icons/logo.png"))
        self.resize(300, 400)

        layout = QVBoxLayout(self)

        # Select "All / None" with inline clickable links
        button_layout = QHBoxLayout()
        select_label = QLabel("Select")
        button_layout.addWidget(select_label)

        # Add clickable links for "All" and "None"
        all_none_label = QLabel('<a href="all">All</a> / <a href="none">None</a>')
        all_none_label.setTextFormat(Qt.RichText)
        all_none_label.setTextInteractionFlags(Qt.TextBrowserInteraction)
        all_none_label.linkActivated.connect(self.handle_link_click)
        button_layout.addWidget(all_none_label)

        # Add the layout to the main layout
        layout.addLayout(button_layout)

        # Add checkboxes for each unique type
        self.checkboxes = []
        for edge_type in unique_types:
            checkbox = QCheckBox(edge_type)
            checkbox.setChecked(True)
            self.checkboxes.append(checkbox)
            layout.addWidget(checkbox)

        # Add OK and Cancel buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def handle_link_click(self, link):
        if link == "all":
            self.select_all()
        elif link == "none":
            self.select_none()

    def select_all(self):
        """Select all types."""
        for checkbox in self.checkboxes:
            checkbox.setChecked(True)

    def select_none(self):
        """Deselect all types."""
        for checkbox in self.checkboxes:
            checkbox.setChecked(False)

    def get_selected_types(self):
        """Return the list of selected types."""
        return [checkbox.text() for checkbox in self.checkboxes if checkbox.isChecked()]