"""
app.py
The main program for booting up the program
"""

# Needed in order to access command line arguments
import sys

from PyQt6.QtWidgets import QApplication

from main_window import MainWindow


app = QApplication(sys.argv)

# Create our main window
window = MainWindow()
window.show()

# Start the program
app.exec()
