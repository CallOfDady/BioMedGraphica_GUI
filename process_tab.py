# process_tab.py

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from data_process import *  # Import all process functions

class ProcessRow(QWidget):
    def __init__(self, entity_type, id_type, file_path, selected_column, parent=None):
        super().__init__(parent)

        self.entity_type = entity_type
        self.id_type = id_type
        self.file_path = file_path
        self.selected_column = selected_column

        # Layout for the row
        self.layout = QVBoxLayout(self)

        # Display entity type and file path
        entity_layout = QHBoxLayout()
        entity_label = QLabel(entity_type)
        self.process_button = QPushButton("Process")
        self.process_button.setFixedWidth(100)

        # Add widgets to the layout
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
            try:
                process_func(self.entity_type, self.id_type, self.file_path, self.selected_column)
                QMessageBox.information(self, "Success", f"Processing {self.entity_type} completed successfully.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to process {self.entity_type}: {str(e)}")
        else:
            QMessageBox.critical(self, "Error", f"No processing function found for entity type: {self.entity_type}")

class ProcessTab(QWidget):
    def __init__(self, read_info_list, parent=None):
        super().__init__(parent)

        self.layout = QVBoxLayout(self)
        self.process_rows = []

        # Create a ProcessRow for each entry in read_info_list
        for entity_type, id_type, file_path, selected_column in read_info_list:
            row = ProcessRow(entity_type, id_type, file_path, selected_column)
            self.process_rows.append(row)
            self.layout.addWidget(row)

        # Add some spacing at the bottom
        self.layout.addStretch()

    def get_process_info(self):
        """Get selected columns from all rows."""
        process_info_list = []
        for row in self.process_rows:
            selected_column = row.selected_column
            process_info_list.append(selected_column)
        return process_info_list
