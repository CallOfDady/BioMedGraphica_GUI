# read_tab.py

import os
import pandas as pd
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

def read_file(file_path, id_type):
    """Read a file and return a DataFrame with appropriate columns based on id_type."""
    _, file_extension = os.path.splitext(file_path)
    
    if file_extension == ".csv":
        df = pd.read_csv(file_path)
    elif file_extension in [".txt", ".tsv"]:
        df = pd.read_csv(file_path, delimiter='\t')
    else:
        raise ValueError(f"Unsupported file type: {file_extension}")
    
    # Check if id_type is "Locus-based ID" to apply different column filtering
    if id_type == "Locus-based ID":
        # Get columns containing "start" or "end"
        start_columns = [col for col in df.columns if "start" in col.lower()]
        end_columns = [col for col in df.columns if "end" in col.lower()]
        
        # Create combinations of start and end columns
        valid_columns = [f"{start}, {end}" for start in start_columns for end in end_columns]
    else:
        # Default behavior: filter columns with "id" in the name
        valid_columns = [col for col in df.columns if "id" in col.lower() or "name" in col.lower()]
        # or col.startswith("Unnamed")
    
    return valid_columns

class ReadRow(QWidget):
    def __init__(self, feature_label, entity_type, id_type, file_path, parent=None):
        super().__init__(parent)

        # Layout for the row
        self.layout = QVBoxLayout(self)
        
        # Display feature label, entity type, and file path
        entity_layout = QHBoxLayout()
        feature_label_display = QLabel(f"Label: {feature_label}")
        entity_label = QLabel(f"Entity Type: {entity_type}")
        self.column_select = QComboBox()
        self.column_select.setFixedWidth(300)  # Adjusted width for the combined columns
        
        # Load file and populate the column dropdown
        try:
            columns = read_file(file_path, id_type)
            self.column_select.addItems(columns)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to read file: {e}")
        
        entity_layout.addWidget(feature_label_display)
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
        for feature_label, entity_type, id_type, file_path in file_info_list:
            row = ReadRow(feature_label, entity_type, id_type, file_path)  # Pass id_type to ReadRow
            self.read_rows.append(row)
            self.layout.addWidget(row)

        # Add some spacing at the bottom
        self.layout.addStretch()

    def get_read_info(self):
        """Get selected columns from all rows, including file_info_list content."""
        read_info_list = []
        for (feature_label, entity_type, id_type, file_path), row in zip(self.file_info_list, self.read_rows):
            selected_column = row.get_selected_column()
            # Combine file_info with the selected column
            read_info_list.append((feature_label, entity_type, id_type, file_path, selected_column))
        return read_info_list