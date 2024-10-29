import os
import shutil
import numpy as np
import traceback
import sys
from itertools import cycle
from sklearn.model_selection import KFold
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

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

class ExportTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.layout = QVBoxLayout(self)

        # Create tree widget for cache folder structure
        self.cache_tree = QTreeWidget()
        self.cache_tree.setHeaderLabel("Cache Folder Structure")
        self.layout.addWidget(self.cache_tree)

        # Populate cache tree
        self.populate_cache_tree()

        # K-Fold Split section
        kfold_label_layout = QVBoxLayout()
        kfold_label_layout.addWidget(QLabel("K-Fold Split:"))
        self.layout.addLayout(kfold_label_layout)

        # K-Fold options layout
        kfold_layout = QHBoxLayout()
        self.kfold_checkbox = QCheckBox()
        self.kfold_combobox = QComboBox()
        self.kfold_combobox.addItems(["3", "5", "10"])
        kfold_layout.addWidget(self.kfold_checkbox)
        kfold_layout.addWidget(self.kfold_combobox)
        kfold_layout.addWidget(QLabel("-Fold Split"))
        self.layout.addLayout(kfold_layout)

        # Export path and controls
        export_layout = QVBoxLayout()
        export_path_label = QLabel("Export Path:")
        export_layout.addWidget(export_path_label)

        # Path input and browse button
        path_layout = QHBoxLayout()
        self.export_path_input = QLineEdit("./output")
        path_layout.addWidget(self.export_path_input)
        self.browse_button = QPushButton("Browse")
        self.browse_button.clicked.connect(self.browse_export_path)
        path_layout.addWidget(self.browse_button)
        export_layout.addLayout(path_layout)

        # Export button
        self.export_button = QPushButton("Export")
        self.export_button.clicked.connect(self.start_export_process)
        export_layout.addWidget(self.export_button)

        # Loading animation
        self.loading_label = QLabel("")
        export_layout.addWidget(self.loading_label)

        self.layout.addLayout(export_layout)

        self.threadpool = QThreadPool()

    def populate_cache_tree(self):
        """Populate tree widget with contents of ./cache directory, skipping folders with '_' prefix."""
        cache_dir = "./cache"
        self.cache_tree.clear()

        if not os.path.exists(cache_dir):
            QMessageBox.warning(self, "Warning", "Cache directory does not exist.")
            return

        for item in os.listdir(cache_dir):
            if item.startswith("_"):
                continue  # Skip folders with '_' prefix
            item_path = os.path.join(cache_dir, item)
            if os.path.isdir(item_path):
                folder_item = QTreeWidgetItem(self.cache_tree, [item])
                self.add_files_to_tree(folder_item, item_path)

        self.cache_tree.expandAll()

    def add_files_to_tree(self, parent_item, folder_path):
        """Recursively add files and subfolders to the tree widget."""
        for item in os.listdir(folder_path):
            if item.startswith("_"):
                continue  # Skip folders with '_' prefix
            item_path = os.path.join(folder_path, item)
            if os.path.isdir(item_path):
                subfolder_item = QTreeWidgetItem(parent_item, [item])
                self.add_files_to_tree(subfolder_item, item_path)
            else:
                QTreeWidgetItem(parent_item, [item])

    def browse_export_path(self):
        """Open a dialog to select the export directory."""
        directory = QFileDialog.getExistingDirectory(self, "Select Export Directory", "./")
        if directory:
            self.export_path_input.setText(directory)

    def start_export_process(self):
        """Begin the export process with loading animation and multithreading."""
        worker = Worker(self.export_cache_data)
        worker.signals.finished.connect(self.on_export_complete)
        worker.signals.error.connect(self.on_export_error)
        worker.signals.result.connect(self.on_export_result)

        # Display loading animation
        self.animation_chars = cycle(["⢿", "⣻", "⣽", "⣾", "⣷", "⣯", "⣟", "⡿"])
        self.loading_timer = QTimer()
        self.loading_timer.timeout.connect(self.update_loading_animation)
        self.loading_timer.start(100)

        self.threadpool.start(worker)

    def update_loading_animation(self):
        """Update loading animation in label."""
        char = next(self.animation_chars)
        self.loading_label.setText(f"Exporting... {char}")

    def export_cache_data(self):
        """Export the contents of ./cache subfolders to the selected export path."""
        cache_dir = "./cache"
        export_dir = self.export_path_input.text()

        if not os.path.exists(cache_dir):
            raise Exception("Cache directory does not exist.")

        if not os.path.exists(export_dir):
            os.makedirs(export_dir, exist_ok=True)

        for item in os.listdir(cache_dir):
            if item.startswith("_"):
                continue
            item_path = os.path.join(cache_dir, item)
            if os.path.isdir(item_path):
                shutil.copytree(item_path, os.path.join(export_dir, item))

        # K-Fold split if checkbox is selected
        if self.kfold_checkbox.isChecked():
            k_folds = int(self.kfold_combobox.currentText())
            self.perform_kfold_split(k_folds, export_dir)

    def perform_kfold_split(self, k_folds, export_dir):
        """Perform K-Fold split on xAll.npy and yAll.npy and save in export_dir."""
        processed_data_dir = os.path.join(export_dir, "processed_data")
        os.makedirs(processed_data_dir, exist_ok=True)

        x_all = np.load("./cache/processed_data/xAll.npy")
        y_all = np.load("./cache/processed_data/yAll.npy")
        kf = KFold(n_splits=k_folds, shuffle=True, random_state=42)

        for fold_index, (train_index, test_index) in enumerate(kf.split(x_all)):
            x_train, x_test = x_all[train_index], x_all[test_index]
            y_train, y_test = y_all[train_index], y_all[test_index]
            
            # Save each fold's data
            np.save(os.path.join(processed_data_dir, f"x_train_fold_{fold_index + 1}.npy"), x_train)
            np.save(os.path.join(processed_data_dir, f"x_test_fold_{fold_index + 1}.npy"), x_test)
            np.save(os.path.join(processed_data_dir, f"y_train_fold_{fold_index + 1}.npy"), y_train)
            np.save(os.path.join(processed_data_dir, f"y_test_fold_{fold_index + 1}.npy"), y_test)

    def on_export_complete(self):
        """Stop the loading animation on completion."""
        self.loading_timer.stop()
        self.loading_label.setText("Export completed successfully.")

    def on_export_error(self, error):
        """Handle errors that occur during export."""
        exctype, value, traceback_str = error
        self.loading_timer.stop()
        self.loading_label.setText("")
        QMessageBox.critical(self, "Error", f"An error occurred: {value}\n{traceback_str}")

    def on_export_result(self, result):
        QMessageBox.information(self, "Success", f"Export completed successfully.")
