import math
import os
import random
from math import exp

import numpy as np


if __name__ == "__main__":

    #checking synthesis efficiency
    """count_i = 0
    with open(fr'{os.getcwd()}\dna-fountain\synthesis_file_2.txt') as file:
        for line in file:
            if len(line.split(",")[1].strip()) == 104 :
                count_i += 1

    print(count_i)"""

    path = fr'{os.getcwd()}\dna-fountain'

    e_i = 0.9
    for c in range(10):
        e_i = e_i - (e_i/10)
        print(e_i)



