import os
import math
import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

fSize = 12

numSam = 100

dataDir = '/Users/waynezheng/Downloads/data/'
fileMag = 'ising_mag_size%s_beta%s_%s_%s.dat'

def LoadData(f):
    rawData = np.fromfile(f, dtype=np.float)
    #  eigVal = np.reshape(rawData, (numSam, numEval))
    mag = []
    magErr = []
    for i in range(numSam):
        j = i*2
        k = i*2+1
        mag = np.append(mag, rawData[j])
        magErr = np.append(magErr, rawData[k])
    return mag, magErr


fig = plt.figure(figsize = (8, 4))
mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('text', usetex=True)
# plt.rc('font', family= 'serif')
x = np.linspace(0.01, 1.0, num=numSam)
print(x)
ax = fig.add_subplot(111)
# ax.set_title('mag', fontsize=fsize, x=0.2, y=0.9)
ax.set_xlim(0.0, 1.0)
ax.set_ylim(0.0, 1.0)
ax.axvline(x=0.4497, color='red', linestyle='-.', linewidth=1)


size = 1010
paras = (size, 0.01, 0.01, 100)
mag, magErr = LoadData(os.path.join(dataDir, fileMag) % paras)
ax.errorbar(x, mag, magErr, color='red', linestyle='-', linewidth=0.5, fmt='o', markersize=2, capsize=2, label= '%s' % size)

size = 1616
paras = (size, 0.01, 0.01, 100)
mag, magErr = LoadData(os.path.join(dataDir, fileMag) % paras)
ax.errorbar(x, mag, magErr, color='blue', linestyle='-', linewidth=0.5, fmt='o', markersize=2, capsize=2, label= '%s' % size)

size = 3232
paras = (size, 0.01, 0.01, 100)
mag, magErr = LoadData(os.path.join(dataDir, fileMag) % paras)
ax.errorbar(x, mag, magErr, color='black', linestyle='-', linewidth=0.5, fmt='o', markersize=2, capsize=2, label= '%s' % size)

# ax.plot(x[1::2], bEE0[1::2], '->', color='blue', label='bEE', linewidth=0.5)
ax.legend(loc='upper right', frameon=False, fontsize=fSize)
# ax.set_xlabel('$\beta$', fontsize=fSize)
# ax1.set_yticklabels('', visible=False)
# plt.setp(ax1.get_yticklabels(), visible=False)
# plt.ylabel('spin current aplitude', fontsize=fs)
# plt.setp(ax1.get_xticklabels(), fontsize=fs)
# ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax.tick_params(axis='both', labelsize=fSize, direction='in')
# ax1.tick_params('x', labelsize=16)
# ax1.tick_params('y', labelsize=16)

image = 'magnetization.pdf'
# paras_image = (size, bounCon)
fig.tight_layout()
plt.savefig(image, format='PDF')
# plt.show()
