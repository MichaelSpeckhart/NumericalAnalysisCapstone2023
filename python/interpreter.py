"""
interpreter.py
Code for interpreting user commands
"""

import re
from cmd_list import cmd_list


def parse_command(args):
    # display matrices
    regex = re.match(r"print\((\d*)\)", args)
    if regex:
        global matrices
        if regex.group(1) == "":
            return [cmd_list.ERROR, "invalid syntax"]
        i = int(regex.group(1))
        return [cmd_list.PRINT, i]

    # sum matrices
    regex = re.match(r"sum\((\d+),(\d+),(y|n)\)", args)
    if regex:
        if regex.group(1) == "" or regex.group(2) == "":
            return [cmd_list.ERROR, "empty matrix index"]
        m1 = int(regex.group(1))
        m2 = int(regex.group(2))
        save = regex.group(3)

        if save == "y":
            return [cmd_list.SUM, m1, m2, True]
        elif save == "n":
            return [cmd_list.SUM, m1, m2, False]
        return [cmd_list.ERROR, "save option must either be 'y' or 'n'"]

    # sub matrices
    regex = re.match(r"sub\((\d+),(\d+),(y|n)\)", args)
    if regex:
        if regex.group(1) == "" or regex.group(2) == "":
            return [cmd_list.ERROR, "empty matrix index"]
        m1 = int(regex.group(1))
        m2 = int(regex.group(2))
        save = regex.group(3)

        if save == "y":
            return [cmd_list.SUB, m1, m2, True]
        elif save == "n":
            return [cmd_list.SUB, m1, m2, False]
        return [cmd_list.ERROR, "save option must either be 'y' or 'n'"]

    # mult matrices
    regex = re.match(r"mult\((\d+),(\d+),(y|n)\)", args)
    if regex:
        if regex.group(1) == "" or regex.group(2) == "":
            return [cmd_list.ERROR, "empty matrix index"]
        m1 = int(regex.group(1))
        m2 = int(regex.group(2))
        save = regex.group(3)

        if save == "y":
            return [cmd_list.MULT, m1, m2, True]
        elif save == "n":
            return [cmd_list.MULT, m1, m2, False]
        return [cmd_list.ERROR, "save option must either be 'y' or 'n'"]

    # rand matrix
    regex = re.match(r"rand\((\d+),(\d+),(\-?\d+),(\-?\d+),(y|n)\)", args)
    if regex:
        if regex.group(1) == "" or regex.group(2) == "":
            return False
        d1 = int(regex.group(1))
        d2 = int(regex.group(2))
        min = int(regex.group(3))
        max = int(regex.group(4))
        save = regex.group(5)

        if min >= max:
            return [cmd_list.ERROR, "invalid bounds"]

        if d1 < 1 or d2 < 1:
            return [cmd_list.ERROR, "invalid dimensions"]

        if save == "y":
            return [cmd_list.RAND, d1, d2, min, max, True]
        elif save == "n":
            return [cmd_list.RAND, d1, d2, min, max, False]
        return [cmd_list.ERROR, "save option must either be 'y' or 'n'"]

    # transpose matrix
    regex = re.match(r"transpose\((\d+),(y|n)\)", args)
    if regex:
        i = int(regex.group(1))
        save = regex.group(2)

        if save == "y":
            return [cmd_list.TRANSPOSE, i, True]
        elif save == "n":
            return [cmd_list.TRANSPOSE, i, False]
        return [cmd_list.ERROR, "save option must either be 'y' or 'n'"]

    # drop matrix
    regex = re.match(r"drop\((\d*)\)", args)
    if regex:
        if regex.group(1) == "":
            return [cmd_list.DROP_ALL]
        i = int(regex.group(1))
        return [cmd_list.DROP, i]

    # plot with colorbar
    regex = re.match(r"plotcolor\((\d+)\)", args)
    if regex:
        i = int(regex.group(1))
        return [cmd_list.PLOT_COLOR, i]

    # descriptive stats
    regex = re.match(r"stats\((\d+),(\d*)\)", args)
    if regex:
        i = int(regex.group(1))
        choice = regex.group(2)

        if choice == "":
            return [cmd_list.STATS, 1, i]
        elif choice == "*":
            return [cmd_list.STATS, 2, i]
        return [cmd_list.STATS, 3, i, int(choice)]

    # 2d scatter plot
    regex = re.match(r"scatter2d\((\d+)\)", args)
    if regex:
        i = int(regex.group(1))
        if i is None or i == "":
            return [cmd_list.ERROR, "must input an index"]
        return [cmd_list.SCATTER2D, i]

    # 3d scatter plot
    regex = re.match(r"scatter3d\((\d+)\)", args)
    if regex:
        i = int(regex.group(1))
        if i is None or i == "":
            return [cmd_list.ERROR, "must input an index"]
        return [cmd_list.SCATTER3D, i]

    # solve linear system
    regex = re.match(r"solve\((\d+)\)", args)
    if regex:
        i = int(regex.group(1))
        if i is None or i == "":
            return [cmd_list.ERROR, "must input an index"]
        return [cmd_list.SOLVE, i]

    return [cmd_list.ERROR, "invalid command"]
