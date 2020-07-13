import numpy as np
import matplotlib.pyplot as plt

def write_dict_as_str(filename, mydictionary, security='low'):
    """
    Writes a dict as a text file. Import requires eval so not particularly safe. 
    """
    
    with open(filename, 'w') as f:
        print(mydictionary, file=f)
    
def read_dict_from_str(filename, security='low'):
    """
    Reads text file produced by textfile_2_dict. Use eval so not particularly safe operation.
    """
    
    with open(filename, 'r') as f:
        line = f.readline().strip()
        mydictionary = eval(line)
    
    return mydictionary
