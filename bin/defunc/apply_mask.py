import numpy as np
import sys


mask = np.loadtxt(sys.argv[3], dtype=bool)
np.save(sys.argv[2], np.load(sys.argv[1])[mask, :])
