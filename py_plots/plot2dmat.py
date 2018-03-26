import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt
import scipy.io as sio
import seaborn as sns


mat_contents = sio.loadmat('matlab.mat')

#x = pd.DataFrame(np.random.randn(9, 9))
x = pd.DataFrame(mat_contents['d'])
x = np.around(x, decimals=2)

samples = [10, 20, 30, 40, 60, 80, 100, 160, 200]
blocks = [2, 4, 8, 13, 15, 17, 25, 30, 40]

# Draw a heatmap with the numeric values in each cell
f, ax = plt.subplots(figsize=(9, 9))
h = sns.heatmap(x, annot=True, fmt=".2f", linewidths=.5, ax=ax, yticklabels=samples, xticklabels=blocks)

h.set_xlabel('number of blocks')
h.set_ylabel('number of samples')

ax.set_title('Random numbers from .mat source file')

plt.show()


