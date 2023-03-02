"""
main_window.py
The main window of the program (really the only window other than plots)
"""
from pandas import *
import numpy as np
import matplotlib

matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import *
from itertools import chain

from PyQt6 import *
from PyQt6.QtWidgets import *
from PyQt6.QtCore import *

import functions

from interpreter import parse_command
from cmd_list import cmd_list
from pandas_model import PandasModel
from table_model import TableModel
from error import *


class MainWindow(QMainWindow):
    """
    Goes through all of the program configuration steps
    """

    def __init__(self):
        super().__init__()
        self.init_ui()
        self.setup_connections()
        self._matrices = []
        self._count = 0

    def init_ui(self):
        self.resize(1300, 700)
        self.center()
        self.setWindowTitle("Capstone Project")
        self.init_layout()

    def center(self):
        qr = self.frameGeometry()
        cp = self.screen().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def init_layout(self):
        self._load_button = QPushButton("Load File")
        self._save_button = QPushButton("Save")
        self._run_button = QPushButton("Run")

        layout = QVBoxLayout()
        button_layout = QHBoxLayout()

        button_layout.addWidget(self._load_button)
        button_layout.addWidget(self._save_button)
        button_layout.addWidget(self._run_button)

        button_layout.insertSpacing(2, 1100)

        self._input_box = QLineEdit()
        self._dialog_box = QTextBrowser()
        self._dialog_box.setFixedHeight(50)
        self._table = QTableView()

        layout.addLayout(button_layout)
        layout.addWidget(self._input_box)
        layout.addWidget(self._dialog_box)
        layout.addWidget(self._table)

        self._central_widget = QWidget()
        self._central_widget.setLayout(layout)
        self.setCentralWidget(self._central_widget)

    """
    Sets up signals and slots for the buttons.
    """

    def setup_connections(self):
        self._load_button.clicked.connect(self.load_button_clicked)
        self._save_button.clicked.connect(self.save_button_clicked)
        self._run_button.clicked.connect(self.run_button_clicked)

    """ 
    Bridge for all of the button clicks.
    """

    def load_button_clicked(self):
        self.execute_load()

    def save_button_clicked(self):
        self.execute_save()

    def run_button_clicked(self):
        self.update_output("")
        self.execute_run()

    """
    Loading a file.
    """

    def execute_load(self):
        file_dialog = QFileDialog.getOpenFileName(
            self, "Open File", "./", "All Files (*)"
        )
        file_name = str(file_dialog[0])
        if file_name is not None and file_name != "":
            read_result = functions.read_file(file_name)
            x = self._count
            if read_result == [[]]:
                self.update_output(
                    "error: either the file was empty or there was an error loading the data"
                )
                return
            for matrix in read_result:
                self.save_matrix(matrix)
            if self._count - x > 1:
                self.update_output(
                    "{} matrices successfully loaded".format(self._count - x)
                )
            else:
                self.update_output("1 matrix successfully loaded")

    """
    Saving a file.
    """

    def execute_save(self):
        file_dialog = QFileDialog.getOpenFileName(
            self, "Open File", "./", "All Files (*)"
        )
        file_name = str(file_dialog[0])
        if file_name is not None and len(self._matrices) > 0:
            if functions.save_file(self._matrices, file_name):
                self.update_output("File successfully saved")
                return
            self.update_output("Failed to save file")

    """
    Running a command.
    """

    def execute_run(self):
        cmd = str(self._input_box.text())
        res = parse_command(cmd)

        if len(res) < 0:
            self.update_output("something went wrong...")
            return
        
    
        match res[0]:
            case cmd_list.ERROR:
                """
                ERROR IN COMMAND
                """
                if len(res) == 2:
                    if res[1] is not None or res[1] != "":
                        self.update_output("error: {}".format(res[1]))
                    return
                self.update_output("unknown error")
            case cmd_list.PRINT:
                """
                PRINT OPERATION
                """
                i = res[1]
                if i is not None:
                    if self.check_index(i):
                        self.display_matrix(self._matrices[i])
                        return
                    self.update_output("invalid index")
                    return
                self.update_output("index error")
            case cmd_list.SUM:
                """
                SUM OPERATION
                """
                m1, m2, save = res[1], res[2], res[3]
                if self.check_index(m1) and self.check_index(m2):
                    self.handle_result(
                        functions.sum_matrix(self._matrices[m1], self._matrices[m2]),
                        save,
                    )
                    return
                self.update_output("index error")
            case cmd_list.SUB:
                """
                SUBTRACT OPERATION
                """
                m1, m2, save = res[1], res[2], res[3]
                if self.check_index(m1) and self.check_index(m2):
                    self.handle_result(
                        functions.sub_matrix(self._matrices[m1], self._matrices[m2]),
                        save,
                    )
                    return
                self.update_output("index error")
            case cmd_list.MULT:
                """
                MULTIPLY OPERATION
                """
                m1, m2, save = res[1], res[2], res[3]
                if self.check_index(m1) and self.check_index(m2):
                    self.handle_result(
                        functions.mult_matrix(self._matrices[m1], self._matrices[m2]),
                        save,
                    )
                    return
                self.update_output("index error")
            case cmd_list.RAND:
                self.handle_result(
                    functions.generate_random_matrix(res[1], res[2], res[3], res[4]),
                    res[5],
                )
            case cmd_list.TRANSPOSE:
                """
                TRANSPOSE OPERATION
                """
                i, save = res[1], res[2]
                if self.check_index(i) is False:
                    self.update_output("index error")
                    return
                self.handle_result(functions.transpose(self._matrices[i]), save)

            case cmd_list.DROP:
                """
                DROP OPERATION
                """
                i = res[1]
                if self.check_index(i):
                    self._matrices.pop(i)
                    self._count -= 1
                    self.update_output("matrix {} dropped".format(i))
                    return
                self.update_output("index error")
            case cmd_list.DROP_ALL:
                """
                DROP ALL OPERATION
                """
                self._matrices.clear()
                self._count = 0
            case cmd_list.PLOT_COLOR:
                """
                COLOR PLOT OPERATION
                """
                i = res[1]
                if self.check_index(i):
                    plt.imshow(np.array(self._matrices[i]))
                    plt.colorbar()
                    plt.show()
                    return
                self.update_output("index error")
            case cmd_list.STATS:
                """
                STATS OPERATION
                """
                choice = res[1]
                i = res[2]
                if self.check_index(i) == False:
                    self.update_output("index error")
                    return
                if choice == 1:
                    # get stats on entire matrix
                    df = (
                        DataFrame(list(chain.from_iterable(self._matrices[i])))
                        .describe()
                        .transpose()
                    )
                    self.display_pd_model(df.values.tolist(), df.columns.tolist())
                elif choice == 2:
                    # get stats for each row
                    matrix = self._matrices[i]
                    frames = []
                    columns = []
                    for row in matrix:
                        df = DataFrame(row).describe().transpose()
                        frames.append(df.values.tolist())
                        columns = df.columns.tolist()
                    print(frames)
                    self.display_pd_model_mult(frames, columns)
                elif choice == 3:
                    # get stats specific row
                    row = res[3]
                    if self.out_of_bounds_matrix(i, row, True):
                        self.update_output("index in matrix out of bounds")
                        return
                    df = DataFrame(self._matrices[i][row]).describe().transpose()
                    self.display_pd_model(df.values.tolist(), df.columns.tolist())
            case cmd_list.SCATTER2D:
                """
                2D SCATTER PLOT OPERATION
                """
                i = res[1]
                if self.check_index(i) is False:
                    self.update_output("invalid index")
                    return
                m = self._matrices[i]
                if len(m[0]) != 2:
                    self.update_output(
                        "invalid dimensions: matrix does not have two columns"
                    )
                    return

                columns = list(zip(*m))
                plt.scatter(columns[0], columns[1], c="blue", marker="o")
                plt.show()
            case cmd_list.SCATTER3D:
                """
                3D SCATTER PLOT OPERATION
                """
                i = res[1]
                if self.check_index(i) is False:
                    self.update_output("invalid index")
                    return
                m = self._matrices[i]
                if len(m[0]) != 3:
                    self.update_output(
                        "invalid dimensions: matrix does not have three columns"
                    )
                    return

                columns = list(zip(*m))
                fig = plt.figure()
                ax = fig.add_subplot(111, projection="3d")
                ax.scatter(columns[0], columns[1], columns[2], c="r", marker="o")
                ax.set_xlabel("x axis")
                ax.set_ylabel("y axis")
                ax.set_zlabel("z axis")
                plt.show()
            case cmd_list.SOLVE:
                """
                SOLVE LINEAR SYSTEM
                """
                i = res[1]
                if self.check_index(i) is False:
                    self.update_output("invalid index")
                m = self._matrices[i]
                l = len(m)

                res = functions.gaussian_elimination(m)

                if len(res[0]) == l + 1:
                    self.update_output("inconsistent system")
                    return
                if len(res[0]) == l + 2:
                    self.update_output("infinately many solutions")
                    return
                if res == d_err:
                    self.update_output("dimensions error")
                self.display_matrix(res)

    """
    Checks for out of bounds error.
    """

    def out_of_bounds_matrix(self, m, i, check_row):
        matrix = self._matrices[m]
        if check_row:
            if i > len(matrix):
                return True
            return False
        if i > len(matrix[0]):
            return True
        return False

    """
    Handles a results by displaying error messages on failure or 
    saving and displaying results on success.
    """

    def handle_result(self, res, save):
        if res == a_err:
            self.update_output(a_err_msg)
        elif res == d_err:
            self.update_output(d_err_msg)
        elif save:
            self.save_matrix(res)
        else:
            self.display_matrix(res)

    """
    Saves a matrix to the working set.
    """

    def save_matrix(self, matrix):
        self._matrices.append(matrix)
        self.update_output("new matrix[{}] added".format(self._count))
        self._count += 1

    """
    Displays a matrix (mostly a bridge function).
    """

    def display_matrix(self, matrix):
        self._model = TableModel(matrix)
        self._table.setModel(self._model)

    """
    Display the pd model using a table view.
    see: https://doc.qt.io/qtforpython-5/PySide2/QtWidgets/QTableView.html
    """

    def display_pd_model(self, data, columns):
        self._model = PandasModel(data, columns)
        self._table.setModel(self._model)

    """
    This is a helper function...
    """

    def display_pd_model_mult(self, data, columns):
        self._model = PandasModel(data[0], columns)
        data.pop(0)
        for i in data:
            print(i)
            self._model._data.append(i)
        self._table.setModel(self._model)

    """
    This is a helper function...
    """

    def update_output(self, message):
        self._dialog_box.setPlainText(message)

    """
    Checks if the index of a matrix is within the working set.
    """

    def check_index(self, i):
        if i >= 0 and i < self._count:
            return True
        return False
