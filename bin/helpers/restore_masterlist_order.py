import sys
import numpy as np

np_array = np.load(sys.argv[1])
argsort = np.loadtxt(sys.argv[2], dtype=int)
restored_array = np_array[argsort]
np.save(sys.argv[3], restored_array)