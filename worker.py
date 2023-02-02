from PyQt6.QtCore import QRunnable, pyqtSlot

"""
Worker Thread

Inherits worker thread setup and signals from QRunnable with basic but effective functionality

:param fn: The function to be execute when spinning up a new thread
:type fn: function

:param args: List of arguments to be passed into the function

NB: this class is not used but could be useful for potential updates to improve rendering.
"""


class Worker(QRunnable):
    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        self._fn = fn
        self._args = args
        self._kwargs = kwargs

    @pyqtSlot()
    def run(self):
        self._fn(*self._args, **self._kwargs)
