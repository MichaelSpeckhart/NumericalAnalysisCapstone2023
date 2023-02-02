"""
pandas_model.py
A class contructed when displaying matrices.
"""

import sys
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import Qt
import pandas as pd
from decimal import Decimal


class PandasModel(QtCore.QAbstractTableModel):
    def __init__(self, data, columns):
        super(PandasModel, self).__init__()
        self._data = data
        self._columns = columns

    def data(self, index, role):
        if role == Qt.ItemDataRole.DisplayRole:
            return self._data[index.row()][index.column()]

    def rowCount(self, index):
        return len(self._data)

    def columnCount(self, index):
        return len(self._data[0])

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.ItemDataRole.DisplayRole:
            if orientation == Qt.Orientation.Horizontal:
                return self._columns[section]
