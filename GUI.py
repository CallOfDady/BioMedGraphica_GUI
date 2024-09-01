# GUI.py

import sys
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from import_tab import ImportTab
from read_tab import ReadTab  # Import the ReadTab class
from process_tab import ProcessTab  # Import the ProcessTab class

class MedGraphica(QMainWindow):
    def __init__(self):
        super().__init__()

        # Set main window properties
        self.setWindowTitle("MedGraphica")
        self.setWindowIcon(QIcon("images/UI/logo.png"))
        self.resize(1000, 400)  # Initial window size

        # Center the window on the screen
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
    
        # Create the central widget and layout
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        self.main_layout = QVBoxLayout(self.central_widget)

        # Create the tab widget
        self.tab_widget = QTabWidget(self)
        self.main_layout.addWidget(self.tab_widget)

        # Create tabs
        self.import_tab = ImportTab()
        self.read_tab = None  # Initialize ReadTab
        self.process_tab = None  # Initialize ProcessTab
        self.export_tab = QWidget()

        # Add tabs to the tab widget
        self.tab_widget.addTab(self.import_tab, "Import")
        self.tab_widget.addTab(QWidget(), "Read")
        self.tab_widget.addTab(QWidget(), "Process")
        self.tab_widget.addTab(self.export_tab, "Export")

        # Create layouts for export tab (currently empty)
        self.export_layout = QVBoxLayout(self.export_tab)

        # Create the bottom button area
        self.button_layout = QHBoxLayout()
        self.main_layout.addLayout(self.button_layout)

        # Create the Previous, Next, and Finish buttons
        self.previous_button = QPushButton("Previous")
        self.next_button = QPushButton("Next")
        self.finish_button = QPushButton("Finish")

        # Add buttons to the button layout
        self.button_layout.addWidget(self.previous_button)
        self.button_layout.addWidget(self.next_button)
        self.button_layout.addWidget(self.finish_button)

        # Connect button signals to slots
        self.previous_button.clicked.connect(self.go_to_previous_tab)
        self.next_button.clicked.connect(self.go_to_next_tab)
        self.finish_button.clicked.connect(self.finish_read)

        # Initialize button states
        self.update_buttons()

    def go_to_previous_tab(self):
        """Go to the previous tab"""
        current_index = self.tab_widget.currentIndex()
        if current_index > 0:
            self.tab_widget.setCurrentIndex(current_index - 1)
        self.update_buttons()

    def go_to_next_tab(self):
        """Go to the next tab based on the current tab"""
        current_index = self.tab_widget.currentIndex()
        
        print(f"Current Tab Index: {current_index}")  # Debug info
        
        if current_index == 0:
            if not self.handle_import_next():
                return  # If there was an error, do not proceed
        elif current_index == 1:
            self.handle_read_next()
        elif current_index == 2:
            self.handle_process_next()
        elif current_index == 3:
            self.handle_export_next()
        
        # After handling, move to the next tab if appropriate
        if current_index < self.tab_widget.count() - 1:
            self.tab_widget.setCurrentIndex(current_index + 1)
            print(f"Switched to Tab Index: {self.tab_widget.currentIndex()}")  # Debug info
        self.update_buttons()

    def handle_import_next(self):
        """Handle the Next button in the Import tab"""
        file_info_list, error_message = self.import_tab.get_all_file_info()
        if error_message:
            QMessageBox.critical(self, "Error", error_message)
            return False  # Return False if there was an error
        
        # Initialize ReadTab with the file info list
        self.read_tab = ReadTab(file_info_list)
        self.tab_widget.removeTab(1)
        self.tab_widget.insertTab(1, self.read_tab, "Read")
        
        # Print file info to console
        print(f"Total files uploaded: {len(file_info_list)}")
        print("File Info:", file_info_list)
        return True
    
    def handle_read_next(self):
        """Handle the Next button in the Read tab"""
        print("Reading the data...")
        read_info = self.read_tab.get_read_info()
        print(f"Selected Columns: {read_info}")
        
        # Initialize ProcessTab with the read info list
        self.process_tab = ProcessTab(read_info)
        self.tab_widget.removeTab(2)
        self.tab_widget.insertTab(2, self.process_tab, "Process")

    def handle_process_next(self):
        """Handle the Next button in the Process tab"""
        print("Processing the data...")
        # TODO: Add data process logic

    def handle_export_next(self):
        """Handle the Next button in the Export tab"""
        print("Preparing to export the data...")
        # TODO

    def finish_read(self):
        """Finish the read and close the application"""
        self.close()

    def update_buttons(self):
        """Update the states of the navigation buttons"""
        current_index = self.tab_widget.currentIndex()
        self.previous_button.setEnabled(current_index > 0)
        self.next_button.setVisible(current_index < self.tab_widget.count() - 1)
        self.finish_button.setVisible(current_index == self.tab_widget.count() - 1)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MedGraphica()
    window.show()
    sys.exit(app.exec_())
