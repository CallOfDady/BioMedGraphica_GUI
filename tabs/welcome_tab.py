import os
import pandas as pd
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class WelcomeTab(QWidget):
    def __init__(self):
        super().__init__()
        self.database_path = self.load_database_path()  # Load the initial database path if it exists

        # Layout for the WelcomeTab
        main_layout = QVBoxLayout(self)

        # Path selection layout
        path_layout = QHBoxLayout()
        path_label = QLabel("Select BioMedGraphica Path:")
        self.path_input = QLineEdit()
        self.path_input.setText(self.database_path or "")  # Set path if it exists
        self.browse_button = QPushButton("Browse")
        self.browse_button.clicked.connect(self.browse_folder)
        
        path_layout.addWidget(path_label)
        path_layout.addWidget(self.path_input)
        path_layout.addWidget(self.browse_button)

        # Add path selection layout to the main layout
        main_layout.addLayout(path_layout)

        # Create layout for software info, bar plot, and entity counts
        content_layout = QHBoxLayout()

        # Software information with fixed size and link
        software_info_layout = QVBoxLayout()

        # Info label with word-wrapping and fixed width
        info_text = QLabel("For more detailed guidance, please visit the ")
        info_text.setWordWrap(True)
        info_text.setFixedWidth(180)
        
        link_label = QLabel('<a href="https://github.com/FuhaiLiAiLab/BioMedGraphica/blob/main/README.md">GitHub repository link</a>')
        link_label.setWordWrap(True)
        link_label.setOpenExternalLinks(True)
        
        # Add info text and link label to layout
        software_info_layout.addWidget(info_text)
        software_info_layout.addWidget(link_label)
        software_info_layout.addStretch()

        # Software info widget with fixed width
        software_info_widget = QWidget()
        software_info_widget.setFixedWidth(200)
        software_info_widget.setLayout(software_info_layout)
        
        # Placeholder for bar plot
        self.figure = plt.Figure()
        self.canvas = FigureCanvas(self.figure)
        
        # Entity counts display
        self.entity_count_text = QTextEdit()
        self.entity_count_text.setReadOnly(True)
        
        # Adding each part to content layout
        content_layout.addWidget(software_info_widget, stretch=3)
        content_layout.addWidget(self.canvas, stretch=4)
        content_layout.addWidget(self.entity_count_text, stretch=1)

        # Add content layout to the main layout
        main_layout.addLayout(content_layout)

        # If a valid path is already loaded, display entity data
        if self.database_path and self.validate_directory_structure(self.database_path):
            self.load_entity_data(self.database_path)

    def load_database_path(self):
        """Check for config.py and load database_path if available."""
        config_path = "config.py"
        if not os.path.exists(config_path):
            # Create an empty config.py file with database_path = None
            with open(config_path, "w") as f:
                f.write("database_path = None\n")
            return None
        else:
            # Load database_path from config.py
            config_data = {}
            with open(config_path, "r") as f:
                exec(f.read(), config_data)
            return config_data.get("database_path")

    def browse_folder(self):
        """Open a dialog to select the BioMedGraphica path and check its contents."""
        folder_path = QFileDialog.getExistingDirectory(self, "Select BioMedGraphica Directory")
        if folder_path:
            self.path_input.setText(folder_path)
            self.database_path = folder_path  # Store the selected database path
            # Validate the folder structure
            if self.validate_directory_structure(folder_path):
                self.save_database_path(folder_path)  # Save path to config.py
                self.load_entity_data(folder_path)  # Load and display entity data
            else:
                QMessageBox.critical(self, "Error", "Confirm database integrity. 'interaction' and 'node' folders are required.")

    def validate_directory_structure(self, folder_path):
        """Check if the directory contains 'interaction' and 'node' folders."""
        expected_folders = {"Interaction", "Node"}
        actual_folders = set(os.listdir(folder_path))
        return expected_folders.issubset(actual_folders)

    def save_database_path(self, folder_path):
        """Write the selected database path to config.py."""
        with open("config.py", "w") as f:
            f.write(f"database_path = '{folder_path}'\n")

    def load_entity_data(self, folder_path):
        """Load entity data from the 'node' folder and display it in a bar plot and text display."""
        node_folder = os.path.join(folder_path, "Node")
        entity_counts = {}

        # Traverse each subfolder in the 'node' folder
        for entity_folder in os.listdir(node_folder):
            entity_path = os.path.join(node_folder, entity_folder)
            if os.path.isdir(entity_path):
                # Assume each folder has one CSV file
                csv_files = [f for f in os.listdir(entity_path) if f.endswith('.csv')]
                if csv_files:
                    csv_file_path = os.path.join(entity_path, csv_files[0])
                    # Read CSV and count rows (excluding header)
                    row_count = self.get_csv_row_count(csv_file_path)
                    if row_count is not None:
                        entity_counts[entity_folder] = row_count

        # Visualize and display entity counts if data is loaded
        if entity_counts:
            self.display_entity_counts(entity_counts)

    def get_csv_row_count(self, csv_file_path):
        """Get the row count of a CSV file by reading only the first column."""
        try:
            df = pd.read_csv(csv_file_path, usecols=[0], encoding="utf-8")
            return len(df)  # Row count excluding header
        except UnicodeDecodeError:
            try:
                df = pd.read_csv(csv_file_path, usecols=[0], encoding="ISO-8859-1")
                return len(df)
            except Exception as e:
                QMessageBox.warning(self, "Error", f"Error reading {csv_file_path}: {e}")
                return None

    def display_entity_counts(self, entity_counts):
        """Display entity counts as a bar plot and detailed text information."""
        # Calculate percentages for each entity
        total_count = sum(entity_counts.values())
        entity_percentages = {entity: (count / total_count) * 100 for entity, count in entity_counts.items()}

        # Update bar plot with percentages
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        entities = list(entity_counts.keys())
        percentages = list(entity_percentages.values())
        
        bars = ax.bar(entities, percentages, color="skyblue")
        ax.set_xlabel("Entities")
        ax.set_ylabel("Percentage (%)")
        ax.set_title("Entity Distribution by Percentage")
        ax.set_xticklabels(entities, rotation=45, ha="right")
        
        # Add percentage labels on top of each bar
        for bar, percentage in zip(bars, percentages):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, height, f"{percentage:.1f}%", ha='center', va='bottom')

        self.canvas.draw()

        # Update entity count text without "rows" suffix
        self.entity_count_text.clear()
        for entity, count in entity_counts.items():
            self.entity_count_text.append(f"{entity}: {count}")

    def get_database_path(self):
        """Return the selected database path."""
        return self.database_path
