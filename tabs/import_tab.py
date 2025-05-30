import os
import pandas as pd
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
        self.add_button.setFixedSize(30, 30)
        self.remove_button.setFixedSize(30, 30)

        self.fill0_checkbox = QCheckBox()
        self.fill0_checkbox.setToolTip("Fill 0")

        # Omics Feature Label input
        self.feature_label = QLineEdit()
        self.feature_label.setFixedWidth(300)
        self.feature_label.setPlaceholderText("Omics Feature Label")

        # Data Type dropdown
        self.data_type = QComboBox()
        self.data_type.setFixedWidth(300)
        self.data_type.addItems(["Gene", "Transcript", "Protein", "Promoter", "Drug", "Disease", "Phenotype", "MicroBiome", "Label"])
        self.data_type.setEditable(True)
        self.data_type.lineEdit().setPlaceholderText("Data Type")
        self.data_type.setCurrentIndex(-1)

        # Entity ID type dropdown
        self.id_types_dict = {
            "Gene": ["Ensembl_Gene_ID", "Locus-based ID", "HGNC_Symbol", "Ensembl_Gene_ID_Version", "HGNC_ID", "OMIM_ID", "NCBI_ID", "RefSeq_ID", "GO_ID"],
            "Transcript": ["Ensembl_Transcript_ID", "Ensembl_Transcript_ID_Version", "Ensembl_Gene_ID", "Reactome_ID", "RefSeq_ID", "RNACentral_ID"],
            "Protein": ["Ensembl_Protein_ID", "Ensembl_Protein_ID_Version", "RefSeq_ID", "Uniprot_ID", "HGNC_Symbol"],
            "Promoter": ["Ensembl_Gene_ID", "HGNC_Symbol", "Ensembl_Gene_ID_Version", "HGNC_ID", "OMIM_ID", "NCBI_ID", "RefSeq_ID", "GO_ID"],
            "Drug": ["PubChem_CID_ID", "PubChem_SID_ID", "CAS_ID", "NDC_ID", "UNII_ID", "InChI_ID", "ChEBI_ID", "DrugBank_ID"],
            "Disease": ["OMIM_ID", "ICD11_ID", "ICD10_ID", "DO_ID", "SnomedCT_ID", "UMLS_ID", "MeSHID", "Mondo_ID"],
            "Phenotype": ["Phenotype_Name", "HPO_ID", "OMIM_ID", "Orpha_ID", "UMLS_ID"], 
            "MicroBiome": ["NCBI_ID", "SILVA_ID", "Greengenes_ID", "RDP_ID", "RNACentral_ID", "GTDB_ID"],
            "Label": [""]
        }

        self.id_type = QComboBox()
        self.id_type.setFixedWidth(400)
        self.id_type.setEditable(True)
        self.id_type.lineEdit().setPlaceholderText("ID Type")
        self.id_type.setEnabled(False)

        # Path display
        self.path_display = QLineEdit()
        self.path_display.setReadOnly(True)
        self.path_display.setPlaceholderText("File Path")

        # Upload and clear buttons
        self.upload_button = QPushButton()
        self.upload_button.setFixedSize(30, 30)
        self.upload_button.setIcon(QIcon("assets/icons/upload.png"))
        self.clear_button = QPushButton()
        self.clear_button.setFixedSize(30, 30)
        self.clear_button.setIcon(QIcon("assets/icons/clear.png"))
    
        # Add widgets to the layout
        self.layout.addWidget(self.add_button)
        self.layout.addWidget(self.remove_button)
        self.layout.addWidget(self.fill0_checkbox)
        self.layout.addWidget(self.feature_label)
        self.layout.addWidget(self.data_type)
        self.layout.addWidget(self.id_type)
        self.layout.addWidget(self.path_display)
        self.layout.addWidget(self.upload_button)
        self.layout.addWidget(self.clear_button)

        # Set callbacks
        self.data_type.currentTextChanged.connect(self.update_id_type)
        self.remove_callback = remove_callback
        self.insert_callback = insert_callback
        self.remove_button.clicked.connect(self.remove_row)
        self.add_button.clicked.connect(self.insert_row)
        self.upload_button.clicked.connect(self.upload_file)
        self.clear_button.clicked.connect(self.clear_path)

    def update_id_type(self):
        """Update the id_type dropdown based on the selected data_type"""
        entity = self.data_type.currentText()
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
        default_dir = "./input_data"
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Select a text-based file", 
            default_dir,  
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
        """Get fill0, feature label, Data Type, id type, file path"""
        return (self.fill0_checkbox.isChecked(),
                self.feature_label.text(), 
                self.data_type.currentText(), 
                self.id_type.currentText(), 
                self.path_display.text())

    def set_file_info(self, fill0, feature_label, data_type, id_type, file_path):
        """Set fill0, feature label, Data Type, id type, file path"""
        self.fill0_checkbox.setChecked(bool(fill0))
        self.feature_label.setText(feature_label)
        self.data_type.setCurrentText(data_type)
        self.id_type.setCurrentText(id_type)
        self.path_display.setText(file_path)



class ImportTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        # Main layout
        self.main_layout = QVBoxLayout(self)
        
        # Track if "Label" has been selected
        self.label_selected = False

        # Config file label and buttons at the top
        self.setup_config_controls()
        self.setup_table_header()

        # Upload rows
        self.upload_rows = []
        self.upload_rows_layout = QVBoxLayout()

        # Scroll area for upload rows
        self.scroll_area = QScrollArea()
        self.scroll_widget = QWidget()
        self.scroll_widget.setLayout(self.upload_rows_layout)
        self.scroll_area.setWidget(self.scroll_widget)
        self.scroll_area.setWidgetResizable(True)

        # Add the scroll area to the main layout
        self.main_layout.addWidget(self.scroll_area)

        # Initial N upload rows
        for _ in range(4):
            self.add_row()

    def setup_table_header(self):
        header_layout = QHBoxLayout()


        add_remove_fill0_label = QLabel("Add/Remove   Fill 0")
        add_remove_fill0_label.setAlignment(Qt.AlignCenter)
        add_remove_fill0_label.setFixedWidth(120)

        feature_label = QLabel("Omics Feature Label")
        feature_label.setFixedWidth(300)

        data_type_label = QLabel("Data Type")
        data_type_label.setFixedWidth(300)

        id_type_label = QLabel("ID Type")
        id_type_label.setFixedWidth(400)

        path_label = QLabel("File Path")
        path_label.setFixedWidth(200)

        file_ops_label = QLabel("Upload/Clear")
        file_ops_label.setFixedWidth(70)

        # Font and color settings for labels
        for label in [add_remove_fill0_label, feature_label, data_type_label, id_type_label, path_label, file_ops_label]:
            label.setStyleSheet("color: gray; font-size: 12px;")

        header_layout.addWidget(add_remove_fill0_label)
        header_layout.addWidget(feature_label)
        header_layout.addWidget(data_type_label)
        header_layout.addWidget(id_type_label)
        header_layout.addWidget(path_label)
        header_layout.addWidget(file_ops_label)

        self.main_layout.addLayout(header_layout)

    def setup_config_controls(self):
        """Setup config file label and buttons in a single row."""
        config_label_layout = QHBoxLayout()
        config_label = QLabel("Config File:")

        small_font = QFont()
        small_font.setPointSize(12)
        config_label.setFont(small_font)
        config_label.setFixedSize(150, 40)

        # Add label to layout
        config_label_layout.addWidget(config_label)

        # Config buttons in the same layout as label
        self.import_button = QPushButton("Import")
        self.import_button.setFixedSize(100, 40)
        self.import_button.setFont(small_font)
        self.import_button.clicked.connect(self.import_config_file)

        self.export_button = QPushButton("Export")
        self.export_button.setFixedSize(100, 40)
        self.export_button.setFont(small_font)
        self.export_button.clicked.connect(self.export_config_file)

        config_label_layout.addWidget(self.import_button)
        config_label_layout.addWidget(self.export_button)
        config_label_layout.addStretch()

        # Add the single row layout to main layout
        self.main_layout.addLayout(config_label_layout)

    def add_row(self, after_row=None):
        """Add a new row to the upload area"""
        row = UploadRow(self, remove_callback=self.remove_row, insert_callback=self.add_row)
        row.data_type.currentTextChanged.connect(self.check_label_selection)  # Connect signal to check Label selection
        if after_row:
            index = self.upload_rows.index(after_row) + 1
            self.upload_rows.insert(index, row)
            self.upload_rows_layout.insertWidget(index, row)
        else:
            self.upload_rows.append(row)
            self.upload_rows_layout.addWidget(row)

    def check_label_selection(self):
        """Ensure only one row can have 'Label' selected as data type."""
        label_rows = [row for row in self.upload_rows if row.data_type.currentText() == "Label"]
        if len(label_rows) > 1:
            # More than one Label selected, revert the latest selection and show a message
            QMessageBox.warning(self, "Selection Error", "Only one 'Label' data type can be selected.")
            for row in label_rows[1:]:  # Revert all except the first Label selection
                row.data_type.setCurrentIndex(-1)

    def remove_row(self, row):
        """Remove a row from the upload area"""
        if len(self.upload_rows) > 1:
            self.upload_rows.remove(row)
            self.upload_rows_layout.removeWidget(row)
            row.deleteLater()

    def get_all_file_info(self):
        """Get file info from all rows, allowing feature_label and id_type to be empty if data_type is Label."""
        file_info_list = []
        for row in self.upload_rows:
            fill0, feature_label, data_type, id_type, file_path = row.get_file_info()
            
            # If data_type is not Label, ensure feature_label and id_type are filled
            if data_type != "Label":
                if not feature_label:
                    return None, "Missing Omics Feature Label"
                if not id_type:
                    return None, "Invalid ID Type"
            
            # Append the row's data to file_info_list
            file_info_list.append((fill0, feature_label, data_type, id_type, file_path))

        return file_info_list, None

    def set_all_file_info(self, file_info_list):
        for row in self.upload_rows:
            self.upload_rows_layout.removeWidget(row)
            row.deleteLater()
        self.upload_rows.clear()

        for file_info in file_info_list:
            fill0, feature_label, data_type, id_type, file_path = file_info
            self.add_row()
            self.upload_rows[-1].set_file_info(fill0, feature_label, data_type, id_type, file_path)

    def import_config_file(self):
        """Import a CSV config file and populate the fields, converting all values to strings."""
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Import Config File", 
            "", 
            "CSV Files (*.csv);;All Files (*)", 
            options=options
        )
        if file_path:
            try:
                df = pd.read_csv(file_path, dtype=str)  # Read all values as strings
                required_columns = ['Fill 0', 'Feature Label', 'Data Type', 'ID Type', 'File Path']
                if set(required_columns).issubset(df.columns):
                    df = df[required_columns].fillna('')
                    df['Fill 0'] = df['Fill 0'].map({'True': True, 'False': False}).fillna(False)
                    file_info_list = df[required_columns].fillna('').values.tolist()  # Fill NaN with empty strings
                    self.set_all_file_info(file_info_list)
                else:
                    print("CSV file does not have the required columns.")
            except Exception as e:
                print(f"Failed to import config file: {e}")

    def export_config_file(self):
        """Export the current form to a CSV config file, converting all values to strings."""
        file_info_list, error = self.get_all_file_info()
        if error:
            print(f"Error: {error}")
            return

        # Convert all entries to strings to ensure consistent export
        file_info_list = [[str(item) for item in row] for row in file_info_list]

        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getSaveFileName(
            self, 
            "Export Config File", 
            "Config.csv", 
            "CSV Files (*.csv);;All Files (*)", 
            options=options
        )
        if file_path:
            try:
                df = pd.DataFrame(file_info_list, columns=['Fill 0', 'Feature Label', 'Data Type', 'ID Type', 'File Path'])
                df.to_csv(file_path, index=False)
            except Exception as e:
                print(f"Failed to export config file: {e}")
