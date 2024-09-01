# read_tab.py

import os
import pandas as pd
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

def read_file(file_path):
    """Read a file and return a DataFrame with appropriate columns."""
    _, file_extension = os.path.splitext(file_path)
    
    if file_extension == ".csv":
        df = pd.read_csv(file_path)
    elif file_extension in [".txt", ".tsv"]:
        df = pd.read_csv(file_path, delimiter='\t')
    else:
        raise ValueError(f"Unsupported file type: {file_extension}")
    
    # Filter columns: those with letters and those whose name contains "id" (case-insensitive)
    valid_columns = [col for col in df.columns if "id" in col.lower()]
    
    return df[valid_columns]

class ReadRow(QWidget):
    def __init__(self, entity_type, file_path, parent=None):
        super().__init__(parent)

        # Layout for the row
        self.layout = QVBoxLayout(self)
        
        # Display entity type and file path
        entity_layout = QHBoxLayout()
        entity_label = QLabel(entity_type)
        self.column_select = QComboBox()
        self.column_select.setFixedWidth(220)
        
        # Load file and populate the column dropdown
        try:
            df = read_file(file_path)
            columns = df.columns
            self.column_select.addItems(columns)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to read file: {e}")
        
        entity_layout.addWidget(entity_label)
        entity_layout.addWidget(self.column_select)
        
        # File path display
        path_label = QLabel(file_path)
        path_label.setStyleSheet("color: gray; font-size: 12px;")
        
        # Add widgets to the layout
        self.layout.addLayout(entity_layout)
        self.layout.addWidget(path_label)

    def get_selected_column(self):
        """Get the selected column from the dropdown."""
        return self.column_select.currentText()

class ReadTab(QWidget):
    def __init__(self, file_info_list, parent=None):
        super().__init__(parent)

        self.file_info_list = file_info_list  # Store the file_info_list for later use
        self.layout = QVBoxLayout(self)
        self.read_rows = []

        # Create a ReadRow for each entry in file_info_list
        for entity_type, id_type, file_path in file_info_list:
            row = ReadRow(entity_type, file_path)
            self.read_rows.append(row)
            self.layout.addWidget(row)

        # Add some spacing at the bottom
        self.layout.addStretch()

    def get_read_info(self):
        """Get selected columns from all rows, including file_info_list content."""
        read_info_list = []
        for (entity_type, id_type, file_path), row in zip(self.file_info_list, self.read_rows):
            selected_column = row.get_selected_column()
            # Combine file_info with the selected column
            read_info_list.append((entity_type, id_type, file_path, selected_column))
        return read_info_list
