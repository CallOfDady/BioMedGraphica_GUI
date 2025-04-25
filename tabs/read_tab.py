import os
import pandas as pd
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from itertools import cycle
from functools import partial
import sys
import traceback
import shutil

def read_file(file_path, id_type):
    """Read a file and return a DataFrame with appropriate columns based on id_type."""
    if not file_path or not os.path.isfile(file_path):
        raise ValueError("Invalid or missing file path.")

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
    valid_columns = [col for col in target_columns if "id" in col.lower() or "name" in col.lower() or "sample" in col.lower() or "patient" in col.lower()]
    
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
        self.column_select.clear()
        self.column_select.addItems(columns)
        if columns == ["(Virtual Node)"]:
            self.column_select.setEnabled(False)
        else:
            self.column_select.setEnabled(True)

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
        self.files_loaded = 0  # Counter for completed loads
        self.total_expected_loads = 0
        self.threadpool = QThreadPool()

        self.loading_task = "Loading columns"

        # Add loading label
        self.loading_label = QLabel("Loading columns...")
        self.loading_label.setStyleSheet("font-weight: bold; color: black;")
        self.loading_label.setAlignment(Qt.AlignLeft | Qt.AlignBottom)

        # Loading animation
        self.animation_chars = cycle(["⢿", "⣻", "⣽", "⣾", "⣷", "⣯", "⣟", "⡿"])
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_loading_animation)
        self.timer.start(100)


        # Create a ReadRow for each entry in file_info_list without columns initially
        for index, (fill0, feature_label, entity_type, id_type, file_path) in enumerate(self.file_info_list):
            row = ReadRow(feature_label, entity_type, id_type, file_path)
            self.read_rows.append(row)
            self.layout.addWidget(row)

            self.total_expected_loads += 1

            print(f"[Debug] Starting worker for index={index}")

            def start_worker(i=index, f0=fill0, fp=file_path, idt=id_type):
                if not f0 and fp:
                    worker = Worker(read_file, fp, idt)
                else:
                    worker = Worker(lambda: ["(Virtual Node)"])

                worker.signals.result.connect(lambda columns, ix=i: self.update_row_columns(ix, columns))
                self.threadpool.start(worker)

            start_worker()
                
                
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

        self.layout.addWidget(self.loading_label, alignment=Qt.AlignLeft | Qt.AlignBottom)

    def update_loading_animation(self):
        """Update the loading label with the next animation character."""
        char = next(self.animation_chars)
        self.loading_label.setText(f"{self.loading_task}... {char}")

    def update_row_columns(self, index, columns):
        if 0 <= index < len(self.read_rows):
            self.read_rows[index].update_columns(columns)

        if not hasattr(self, "_loaded_indices"):
            self._loaded_indices = set()

        if index in self._loaded_indices:
            print(f"[Debug] Skipping duplicate update for index={index}")
            return

        self._loaded_indices.add(index)
        self.files_loaded += 1
        print(f"[Debug] update_row_columns called for index={index}, columns={columns}, loaded={self.files_loaded}/{self.total_expected_loads}")

        if self.files_loaded == self.total_expected_loads:
            self.timer.stop()
            self.loading_label.setText("Columns loaded.")

    def get_read_info(self):
        """Return the full file_info + selected_column from each row."""
        read_info_list = []
        for file_info, row in zip(self.file_info_list, self.read_rows):
            selected_column = row.get_selected_column()
            if len(file_info) == 6:
                file_info = (*file_info, selected_column)
            elif len(file_info) == 5:
                file_info = (*file_info, selected_column, list(getattr(self, "intersection_samples", [])))
            read_info_list.append(file_info)
        return read_info_list

    def start_compute_intersection(self):
        """Start the intersection computation in a separate thread."""
        self.intersection_button.setEnabled(False)
        self.loading_task = "Computing intersection"  # Set task type
        self.loading_label.setText(f"{self.loading_task}...")
        self.timer.start(100)

        # Clear cache folders starting with "_"
        self.clear_cache_folders()

        # Create a worker for computing the intersection
        worker = Worker(self.compute_intersection)
        print(f"[Debug] Starting intersection worker")
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
        self.loading_task = "Idle"  # Reset task type
        self.common_samples_count.setText(str(num_common_samples))

    def on_compute_error(self, error):
        """Handle errors during the intersection computation."""
        self.timer.stop()
        self.loading_label.setText("Error during computation.")
        self.loading_task = "Idle"
        self.intersection_button.setEnabled(True)

        exctype, value, tb = error
        print("[ERROR]", exctype, value)
        print(tb)

        QMessageBox.critical(self, "Error", f"{value}\n\n{tb}")

    def clear_cache_folder(self):
        """Clear the cache folder and recreate required subdirectories."""
        cache_folder = './cache'
        raw_id_mapping_dir = os.path.join(cache_folder, "raw_id_mapping")
        processed_data_dir = os.path.join(cache_folder, "processed_data")

        # Delete cache folder and its contents
        if os.path.exists(cache_folder):
            try:
                shutil.rmtree(cache_folder)  # Remove entire cache folder
            except Exception as e:
                print(f"Failed to clear cache folder: {e}")

        # Recreate the cache and required subdirectories
        os.makedirs(raw_id_mapping_dir, exist_ok=True)
        os.makedirs(processed_data_dir, exist_ok=True)

    def compute_intersection(self):
        """Compute the intersection of the first column across all files and save new files."""
        # Clear and recreate cache folder structure
        self.clear_cache_folder()
        
        intersection_samples = None  # For accumulating intersection

        # Load each file's first column and compute intersection
        for file_info in self.file_info_list:
            fill0, feature_label, entity_type, id_type, file_path = file_info
            if fill0:
                continue  # Skip if fill0 is True

            file_path = file_info[4]  # The original file path
            if not os.path.exists(file_path):
                raise ValueError(f"File not found: {file_path}")
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

        # Create output folders for Label and non-Label data
        x_cache_folder = f"./cache/_x_{num_common_samples}"
        y_cache_folder = f"./cache/_y_{num_common_samples}"
        os.makedirs(x_cache_folder, exist_ok=True)
        os.makedirs(y_cache_folder, exist_ok=True)

        # Save new CSV files with intersection samples only
        for file_info in self.file_info_list:
            fill0, feature_label, entity_type, id_type, file_path = file_info

            if fill0 or not file_path.strip():
                continue

            df = pd.read_csv(file_path)

            df = df[df.iloc[:, 0].isin(intersection_samples)]

            if entity_type == "Label":
                df.rename(columns={df.columns[0]: "Sample_ID"}, inplace=True)

            base_name = os.path.basename(file_path)
            new_file_path = os.path.join(
                y_cache_folder if entity_type == "Label" else x_cache_folder,
                f"{base_name.replace('.csv', '')}_{num_common_samples}.csv"
            )

            df.to_csv(new_file_path, index=False)
            file_info[4] = new_file_path  # Update new path
        
        # Add intersection sample list to file_info for downstream usage
        for file_info in self.file_info_list:
            file_info.append(sorted(intersection_samples))

        print(f"[Debug] Intersection samples: {intersection_samples}")

        return num_common_samples