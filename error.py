"""
error.py
All of the error codes and their messages. 
Also defines integer limits on 64 bit machines.
"""

INT_MAX = 2147483647
INT_MIN = -2147483648
d_err = [[INT_MAX, INT_MAX, INT_MIN], [INT_MIN, INT_MIN, INT_MAX]]
a_err = [[INT_MIN, INT_MIN, INT_MAX], [INT_MAX, INT_MAX, INT_MIN]]
d_err_msg = "error in dimensions"
a_err_msg = "arithmetic error"
