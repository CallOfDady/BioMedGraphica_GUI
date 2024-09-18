# process_tab.py

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from itertools import cycle
from data_process import *  # Import all process functions

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
        self.process_button.setFixedWidth(100)

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
    def __init__(self, read_info_list, parent=None):
        super().__init__(parent)

        self.layout = QVBoxLayout(self)
        self.process_rows = []

        # Add loading label
        self.loading_label = QLabel("Processing data...")
        self.layout.addWidget(self.loading_label)
        self.animation_chars = cycle(["⢿", "⣻", "⣽", "⣾", "⣷", "⣯", "⣟", "⡿"])
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_loading_animation)
        self.loading_label.hide()  # Initially hide the label

        # Create a ProcessRow for each entry in read_info_list
        for feature_label, entity_type, id_type, file_path, selected_column in read_info_list:
            row = ProcessRow(feature_label, entity_type, id_type, file_path, selected_column)
            self.process_rows.append(row)
            self.layout.addWidget(row)

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
        worker = Worker(process_func, row.entity_type, row.id_type, row.file_path, row.selected_column, row.feature_label)
        worker.signals.finished.connect(self.on_processing_complete)
        worker.signals.error.connect(self.on_processing_error)
        self.threadpool.start(worker)

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