from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import sys

plt.style.use("classic")

dat = np.loadtxt(sys.argv[1])
Nkpt = dat.shape[0]
Nstates = dat.shape[1] - 1
print("Nkpt = ", Nkpt)
print("Nstates = ", Nstates)

plt.clf()
x = dat[:,0]
for ist in range(1,Nstates+1):
    plt.plot( x, dat[:,ist], marker='o' )

filplot = sys.argv[1].replace(".dat",".pdf")
plt.savefig(filplot)

