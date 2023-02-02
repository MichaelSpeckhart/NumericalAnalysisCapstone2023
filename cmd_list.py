"""
enum.py
Contains an enumeration for all of the commands.
"""

from enum import Enum


class cmd_list(Enum):
    ERROR = 1
    PRINT = 2
    SUM = 3
    SUB = 4
    MULT = 5
    RAND = 6
    TRANSPOSE = 7
    DROP = 8
    DROP_ALL = 9
    PLOT_COLOR = 10
    STATS = 11
    SCATTER2D = 12
    SCATTER3D = 13
    SOLVE = 14
