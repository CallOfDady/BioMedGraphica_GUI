# import_tab.py

import os
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

class UploadRow(QWidget):
    def __init__(self, parent=None, remove_callback=None, insert_callback=None):
        super().__init__(parent)

        # Layout for the row
        self.layout = QHBoxLayout(self)
        
        # Add and remove buttons
        self.add_button = QPushButton("+")
        self.remove_button = QPushButton("-")

        # Set fixed size for add and remove buttons to be square
        self.add_button.setFixedSize(30, 30)
        self.remove_button.setFixedSize(30, 30)

        # Omics Feature Label input
        self.feature_label = QLineEdit()
        self.feature_label.setFixedWidth(150)
        self.feature_label.setPlaceholderText("Omics Feature Label")

        # Entity type dropdown
        self.entity_type = QComboBox()
        self.entity_type.setFixedWidth(150)
        self.entity_type.addItems(["Gene", "Transcript", "Protein", "Drug", "Disease", "Phenotype"])
        self.entity_type.setEditable(True)
        self.entity_type.lineEdit().setPlaceholderText("Entity Type")
        self.entity_type.setCurrentIndex(-1)  # Ensure no item is selected by default

        # Entity id type dropdown
        self.id_types_dict = {
            "Gene": ["Ensembl_Gene_ID", "Locus-based ID", "HGNC_Symbol", "Ensembl_Gene_ID_Version", "HGNC_ID", "OMIM_ID", "NCBI_ID", "RefSeq_ID", "GO_ID"],
            "Transcript": ["Ensembl_Transcript_ID", 'ensembl_gene_id', "Ensembl_Transcript_ID_Version", "Ensembl_Gene_ID", "Reactome_ID", "RefSeq_ID", "RNACentral_ID"],
            "Protein": ["Ensembl_Protein_ID", "Ensembl_Protein_ID_Version", "RefSeq_ID", "Uniprot_ID"],
            "Drug": ["PubChem_CID_ID", "PubChem_SID_ID", "CAS_ID", "NDC_ID", "UNII_ID", "InChI_ID", "ChEBI_ID", "DrugBank_ID"],
            "Disease": ["OMIM_ID", "ICD11_ID", "ICD10_ID", "DO_ID", "SnomedCT_ID", "UMLS_ID", "MeSHID", "Mondo_ID"],
            "Phenotype": ["HPO_ID", "OMIM_ID", "Orpha_ID", "UMLS_ID"], 
            "MicroBiome": ["NCBI_ID", "SILVA_ID", "Greengenes_ID", "RDP_ID", "RNACentral_ID", "GTDB_ID"],
        }

        self.id_type = QComboBox()
        self.id_type.setFixedWidth(220)
        self.id_type.setEditable(True)
        self.id_type.lineEdit().setPlaceholderText("ID Type")
        self.id_type.setEnabled(False)  # Initially disabled

        # Path display
        self.path_display = QLineEdit()
        self.path_display.setReadOnly(True)
        self.path_display.setPlaceholderText("File Path")

        # Upload and clear buttons (also square)
        self.upload_button = QPushButton()
        self.upload_button.setFixedSize(30, 30)
        self.upload_button.setIcon(QIcon("images/UI/upload.png"))

        self.clear_button = QPushButton()
        self.clear_button.setFixedSize(30, 30)
        self.clear_button.setIcon(QIcon("images/UI/clear.png"))
    
        # Add widgets to the layout
        self.layout.addWidget(self.add_button)
        self.layout.addWidget(self.remove_button)
        self.layout.addWidget(self.feature_label)  # Add Omics Feature Label input before Entity Type
        self.layout.addWidget(self.entity_type)
        self.layout.addWidget(self.id_type)
        self.layout.addWidget(self.path_display)
        self.layout.addWidget(self.upload_button)
        self.layout.addWidget(self.clear_button)

        # Set callbacks
        self.entity_type.currentTextChanged.connect(self.update_id_type)
        self.remove_callback = remove_callback
        self.insert_callback = insert_callback
        self.remove_button.clicked.connect(self.remove_row)
        self.add_button.clicked.connect(self.insert_row)
        self.upload_button.clicked.connect(self.upload_file)
        self.clear_button.clicked.connect(self.clear_path)

    def update_id_type(self):
        """Update the id_type dropdown based on the selected entity_type"""
        entity = self.entity_type.currentText()
        self.id_type.clear()

        if entity in self.id_types_dict:
            id_types = self.id_types_dict[entity]
            if id_types:
                self.id_type.addItems(id_types)
                self.id_type.setEnabled(True)
            else:
                self.id_type.setEnabled(False)
        else:
            self.id_type.setEnabled(False)

    def upload_file(self):
        """Open a file dialog to select a file"""
        options = QFileDialog.Options()
        
        # Set a default directory path (you can change this to any directory you prefer)
        default_dir = "E:\LabWork\mosGraphFlow\ROSMAP-raw"  # Default path, modify as needed
        
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Select a text-based file", 
            default_dir,  # This sets the default directory
            "Text-based Files (*.txt *.csv *.tsv);;All Files (*)", 
            options=options
        )
        if file_path:
            self.path_display.setText(file_path)

    def clear_path(self):
        """Clear the file path"""
        self.path_display.clear()

    def remove_row(self):
        """Remove the current row"""
        if self.remove_callback:
            self.remove_callback(self)

    def insert_row(self):
        """Insert a new row below the current one"""
        if self.insert_callback:
            self.insert_callback(self)

    def get_file_info(self):
        """Get feature label, entity type, id type, file path"""
        return (self.feature_label.text(), 
                self.entity_type.currentText(), 
                self.id_type.currentText(), 
                self.path_display.text())

class ImportTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.layout = QVBoxLayout(self)
        self.upload_rows = []

        # Initial 4 upload rows
        for _ in range(2):
            self.add_row()

    def add_row(self, after_row=None):
        """Add a new row to the upload area"""
        row = UploadRow(self, remove_callback=self.remove_row, insert_callback=self.add_row)
        if after_row:
            # Insert the row after the specified row
            index = self.upload_rows.index(after_row) + 1
            self.upload_rows.insert(index, row)
            self.layout.insertWidget(index, row)
        else:
            # Add to the end if no specific row is mentioned
            self.upload_rows.append(row)
            self.layout.addWidget(row)

    def remove_row(self, row):
        """Remove a row from the upload area"""
        if len(self.upload_rows) > 1:
            self.upload_rows.remove(row)
            self.layout.removeWidget(row)
            row.deleteLater()

    def get_all_file_info(self):
        """Get file info from all rows, including feature label"""
        file_info_list = []
        for row in self.upload_rows:
            feature_label, entity_type, id_type, file_path = row.get_file_info()  # Now unpack four values
            if not feature_label:
                return None, "Missing Omics Feature Label"
            if not entity_type:
                return None, "Invalid Entity Type"
            if not id_type:
                return None, "Invalid ID Type"  # Add check for ID Type
            if not file_path or not os.path.exists(file_path):
                return None, "Invalid File Path"
            file_info_list.append((feature_label, entity_type, id_type, file_path))  # Include feature label in the file info
        return file_info_list, None