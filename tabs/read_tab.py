import os
import pandas as pd
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from itertools import cycle
import sys
import traceback
import shutil

def read_file(file_path, id_type):
    """Read a file and return a DataFrame with appropriate columns based on id_type."""
    _, file_extension = os.path.splitext(file_path)
    
    if file_extension == ".csv":
        df = pd.read_csv(file_path, nrows=0)
    elif file_extension in [".txt", ".tsv"]:
        df = pd.read_csv(file_path, delimiter='\t', nrows=0)
    else:
        raise ValueError(f"Unsupported file type: {file_extension}")
        
    first_5_columns = df.columns[:5]
    last_5_columns = df.columns[-5:]
    target_columns = list(first_5_columns) + list(last_5_columns)
    valid_columns = [col for col in target_columns if "id" in col.lower() or "name" in col.lower()]
    
    return valid_columns

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

class ReadRow(QWidget):
    def __init__(self, feature_label, entity_type, id_type, file_path, columns=None, parent=None):
        super().__init__(parent)

        # Layout for the row
        self.layout = QVBoxLayout(self)
        
        # Display feature label, entity type, and file path
        entity_layout = QHBoxLayout()
        feature_label_display = QLabel(f"Label: {feature_label}")
        entity_label = QLabel(f"Entity Type: {entity_type}")
        self.column_select = QComboBox()
        self.column_select.setFixedWidth(300)
        
        # Add widgets to the layout
        entity_layout.addWidget(feature_label_display)
        entity_layout.addWidget(entity_label)
        entity_layout.addWidget(self.column_select)
        
        # File path display
        path_label = QLabel(file_path)
        path_label.setStyleSheet("color: gray; font-size: 12px;")
        
        self.layout.addLayout(entity_layout)
        self.layout.addWidget(path_label)

        # Populate the column dropdown if columns are provided
        if columns:
            self.update_columns(columns)

    def update_columns(self, columns):
        """Update the column dropdown with the provided columns."""
        self.column_select.clear()
        self.column_select.addItems(columns)

    def get_selected_column(self):
        """Get the selected column from the dropdown."""
        return self.column_select.currentText()

class ReadTab(QWidget):
    def __init__(self, file_info_list, parent=None):
        super().__init__(parent)

        # Convert file_info_list to a mutable list
        self.file_info_list = [list(info) for info in file_info_list]
        self.layout = QVBoxLayout(self)
        self.read_rows = []
        self.files_loaded = 0  # Counter for files loaded
        self.threadpool = QThreadPool()

        # Add loading label
        self.loading_label = QLabel("Loading columns...")
        self.layout.addWidget(self.loading_label)
        self.animation_chars = cycle(["⢿", "⣻", "⣽", "⣾", "⣷", "⣯", "⣟", "⡿"])
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_loading_animation)
        self.timer.start(100)

        # Create a ReadRow for each entry in file_info_list without columns initially
        for feature_label, entity_type, id_type, file_path in self.file_info_list:
            row = ReadRow(feature_label, entity_type, id_type, file_path)
            self.read_rows.append(row)
            self.layout.addWidget(row)

        # Add some spacing at the bottom
        self.layout.addStretch()

        # Common samples label and count in the same row
        common_samples_layout = QHBoxLayout()
        self.common_samples_label = QLabel("Common Samples: ")
        self.common_samples_count = QLabel("")  # Initially empty
        common_samples_layout.addWidget(self.common_samples_label)
        common_samples_layout.addWidget(self.common_samples_count)
        self.layout.addLayout(common_samples_layout)

        # Intersection button
        self.intersection_button = QPushButton("Compute Intersection")
        self.intersection_button.clicked.connect(self.start_compute_intersection)
        self.layout.addWidget(self.intersection_button)

    def update_loading_animation(self):
        """Update the loading label with the next animation character."""
        char = next(self.animation_chars)
        self.loading_label.setText(f"Loading... {char}")

    def update_row_columns(self, index, columns):
        """Update the columns for a specific row."""
        if 0 <= index < len(self.read_rows):
            self.read_rows[index].update_columns(columns)
        
        # Increase the count of loaded files
        self.files_loaded += 1

        # Stop the loading animation once all rows are updated
        if self.files_loaded == len(self.read_rows):
            self.timer.stop()
            self.loading_label.setText("Columns loaded.")

    def get_read_info(self):
        """Get selected columns from all rows, including file_info_list content."""
        read_info_list = []
        for (feature_label, entity_type, id_type, file_path), row in zip(self.file_info_list, self.read_rows):
            selected_column = row.get_selected_column()
            read_info_list.append((feature_label, entity_type, id_type, file_path, selected_column))
        return read_info_list

    def start_compute_intersection(self):
        """Start the intersection computation in a separate thread."""
        self.intersection_button.setEnabled(False)
        self.loading_label.setText("Computing intersection...")
        self.timer.start(100)

        # Clear cache folders starting with "_"
        self.clear_cache_folders()

        # Create a worker for computing the intersection
        worker = Worker(self.compute_intersection)
        worker.signals.result.connect(self.on_compute_finished)
        worker.signals.error.connect(self.on_compute_error)
        worker.signals.finished.connect(lambda: self.intersection_button.setEnabled(True))
        self.threadpool.start(worker)

    def clear_cache_folders(self):
        """Remove all folders in the cache directory that start with an underscore."""
        cache_dir = "./cache"
        for folder_name in os.listdir(cache_dir):
            folder_path = os.path.join(cache_dir, folder_name)
            if os.path.isdir(folder_path) and folder_name.startswith("_"):
                try:
                    # Remove the directory
                    shutil.rmtree(folder_path)
                    print(f"Removed cache folder: {folder_path}")
                except Exception as e:
                    print(f"Failed to remove folder {folder_path}: {e}")

    def on_compute_finished(self, num_common_samples):
        """Handle completion of the intersection computation."""
        self.timer.stop()
        self.loading_label.setText("Intersection complete.")
        self.common_samples_count.setText(str(num_common_samples))

    def on_compute_error(self, error):
        """Handle errors during the intersection computation."""
        self.timer.stop()
        exctype, value, traceback_str = error
        QMessageBox.critical(self, "Error", f"An error occurred: {value}\n{traceback_str}")
        self.loading_label.setText("Error during computation.")
        self.intersection_button.setEnabled(True)

    def compute_intersection(self):
        """Compute the intersection of the first column across all files and save new files."""
        intersection_samples = None  # For accumulating intersection

        # Load each file's first column and compute intersection
        for file_info in self.file_info_list:
            file_path = file_info[3]  # The original file path
            df = pd.read_csv(file_path, usecols=[0], dtype=str)
            current_samples = set(df.iloc[:, 0])

            if intersection_samples is None:
                intersection_samples = current_samples
            else:
                intersection_samples &= current_samples

            if not intersection_samples:
                break  # If intersection is empty, no need to continue

        # Count of common samples (excluding header)
        num_common_samples = len(intersection_samples)

        # Create output folder named "./cache/_{num_common_samples}"
        cache_folder = f"./cache/_{num_common_samples}"
        os.makedirs(cache_folder, exist_ok=True)

        # Save new CSV files with intersection samples only
        for file_info in self.file_info_list:
            feature_label, entity_type, id_type, file_path = file_info[:4]
            df = pd.read_csv(file_path)
            df_intersection = df[df.iloc[:, 0].isin(intersection_samples)]
            
            # New filename with row count appended
            base_name = os.path.basename(file_path)
            new_file_path = os.path.join(cache_folder, f"{base_name.replace('.csv', '')}_{num_common_samples}.csv")
            df_intersection.to_csv(new_file_path, index=False)
            
            # Update file path in file_info_list
            file_info[3] = new_file_path  # Update the file path to the new file

        return num_common_samples
