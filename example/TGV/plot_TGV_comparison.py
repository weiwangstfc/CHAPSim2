# this python code is created by W Wang (STFC) to plot data of TGV benchmark
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import rcParams

# Define colormaps and markers
cbrg = cm.get_cmap(name='brg', lut=None)
cgry = cm.get_cmap(name='gray', lut=None)
cbow = cm.get_cmap(name='gist_rainbow', lut=None)
mlst = ["o", "<", "*", "v", "^", '>', '1', '2', '3', '4', 'x', 's', '8', '+']

# Set plot configurations
plt.rc('figure', facecolor="white")
plt.rc('legend', fontsize=15)
rcParams['legend.loc'] = 'best'

figsize = (9, 6)
dpi = 500


ref_path='reference/'
ref_file=ref_path + 'TGV_Re1600.dat'
ref_data=np.genfromtxt(ref_file, delimiter=None, skip_header=43, skip_footer=0,
         names=['t', 'E_k', 'epsilon_t', 'epsilon', 'Dzeta',
                'u2','v2', 'w2', 'dudx2','dudy2',
                'dudz2', 'dvdx2', 'dvdy2', 'dvdz2' , 'dwdx2', 
                'dwdy2', 'dwdz2' ])

dns_path='./'
dns_file='domain1_monitor_bulk.dat'
dns_data=np.genfromtxt(dns_file, delimiter=None, skip_header=2, skip_footer=0,
         names=['t', 'E_k'])

fig1, ax1 = plt.subplots(figsize=figsize, dpi=dpi)
ax1.set_xlabel(r'$time(s)$', fontsize=20)
ax1.set_ylabel(r'$E_k$', fontsize=20)

ax1.plot(ref_data['t'], ref_data['E_k'], marker=mlst[1], mfc='none', ms=4, markevery=10, mec=cbrg(0.00), color=cbrg(0.00), linestyle='None', label=f'ref')
ax1.plot(dns_data['t'], dns_data['E_k'], marker='none', color=cbrg(0.50), linestyle='-.', label=f'CHAPsim2')
ax1.legend(loc='upper right', ncol=1, labelspacing=0.1, frameon=False, handlelength=3.2, numpoints=1)
ax1.grid(color='gray', linestyle='dashed')
fig1.savefig(f'TGV_E_k.png')
print(f"TGV_E_k... Done!")
plt.close('all')












