import numpy as np
import pylab as pl
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.axes_grid1 import host_subplot

plt.rc('figure', facecolor="white")
plt.rc('legend', fontsize=15)
rcParams['legend.loc'] = 'best'
cbrg = cm.get_cmap(name='brg', lut=None) ## 1=green 0=blue  0.5=red
cgry = cm.get_cmap(name='gray',lut=None) ## 1=white 0=black 0.5=grey
#cbow = cm.get_cmap(name='rainbow',lut=None) ## 
cbow = cm.get_cmap(name='gist_rainbow',lut=None) ## 
chsv = cm.get_cmap(name='hsv',lut=None) ## 
csei = cm.get_cmap(name='seismic',lut=None) ## 
ldsh = [(12,0.1),(12,7),(12,5,2,5),(12,5,2,5,2,5),(12,5,2,5,2,5,2,5),(12,7,12,5,2,5), \
        (12,7,12,5,2,5,2,5), (12,7,12,7,12,5,2,5), (12,7,12,7,12,5,2,5,2,5),(12,7),(12,5,2,5),(12,5,2,5,2,5)]
mlst = ["o", "<", "*", "v", "^"]



# load data
Expt=np.genfromtxt('NIST_WATER_23.5MP.DAT', skip_header=2, \
                     skip_footer=0, names=['P', 'H', 'T', 'D', 'M', 'K', 'Cp', 'Beta'])
Sim1=np.genfromtxt('check_tp_from_dh_dim.dat', skip_header=1, \
                     skip_footer=0, names=['H', 'T', 'D', 'M', 'K', 'Cp', 'Beta', 'rhoh'])
Sim2=np.genfromtxt('check_tp_from_dh_undim.dat', skip_header=1, \
                     skip_footer=0, names=['H', 'T', 'D', 'M', 'K', 'Cp', 'Beta', 'rhoh'])

# plot T-H
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Enthalpy$',     fontsize=20)

plt.plot(Expt['T'], Expt['H'],  marker=mlst[0], mfc='none', ms=6,  markevery=10, mec=cbrg(0.0), color=cbrg(0.0), label=r'$exp$')
plt.plot(Sim1['T'], Sim1['H'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_T_H.png')

# plot T-D
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Density$',     fontsize=20)

plt.plot(Expt['T'], Expt['D'],  marker=mlst[0], mfc='none', ms=6,  markevery=10, mec=cbrg(0.0), color=cbrg(0.0), label=r'$exp$')
plt.plot(Sim1['T'], Sim1['D'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_T_D.png')

print("Plotting Check_T_D... Done!")

# plot T-M
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Dynamic Viscosity$',     fontsize=20)

plt.plot(Expt['T'], Expt['M'],  marker=mlst[0], mfc='none', ms=6,  markevery=10, mec=cbrg(0.0), color=cbrg(0.0), label=r'$exp$')
plt.plot(Sim1['T'], Sim1['M'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_T_M.png')

print("Plotting Check_T_M... Done!")

# plot T-K
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Thermal Conductivity$',     fontsize=20)

plt.plot(Expt['T'], Expt['K'],  marker=mlst[0], mfc='none', ms=6,  markevery=10, mec=cbrg(0.0), color=cbrg(0.0), label=r'$exp$')
plt.plot(Sim1['T'], Sim1['K'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_T_K.png')

print("Plotting Check_T_K... Done!")

# plot T-Cp
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Cp$',     fontsize=20)

plt.plot(Expt['T'], Expt['Cp'],  marker=mlst[0], mfc='none', ms=6,  markevery=10, mec=cbrg(0.0), color=cbrg(0.0), label=r'$exp$')
plt.plot(Sim1['T'], Sim1['Cp'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_T_Cp.png')
print("Plotting Check_T_Cp... Done!")

# plot T-rhoh
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Rho*H$',     fontsize=20)

plt.plot(Sim1['T'], Sim1['rhoh'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_T_rhoh.png')

print("Plotting Check_T_rhoh... Done!")


# plot T-H
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Enthalpy$',     fontsize=20)

plt.plot(Sim2['T'], Sim2['H'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_undim_T_H.png')

# plot T-D
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Density$',     fontsize=20)

plt.plot(Sim2['T'], Sim2['D'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_undim_T_D.png')

print("Plotting Check_undim_T_D... Done!")

# plot T-M
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Dynamic Viscosity$',     fontsize=20)

plt.plot(Sim2['T'], Sim2['M'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_undim_T_M.png')

print("Plotting Check_undim_T_M... Done!")

# plot T-K
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Thermal Conductivity$',     fontsize=20)

plt.plot(Sim2['T'], Sim2['K'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_undim_T_K.png')

print("Plotting Check_undim_T_K... Done!")

# plot T-Cp
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Cp$',     fontsize=20)

plt.plot(Sim2['T'], Sim2['Cp'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_undim_T_Cp.png')
print("Plotting Check_undim_T_Cp... Done!")

# plot T-rhoh
fig1 =plt.figure(dpi=500,figsize=(7.5,5.4)) 

pl.xlabel(r'$Temperature $', fontsize=20)
pl.ylabel(r'$Rho*H$',     fontsize=20)

plt.plot(Sim2['T'], Sim2['rhoh'],  marker=mlst[1], mfc='none', ms=6,  markevery=10, mec=cbrg(1.0), color=cbrg(1.0), label=r'$sim$')

plt.legend(loc='upper right',frameon=False,handlelength=3.2, numpoints=1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#gca().ticklabel_format(style='sci',axis='y',scilimits=(-2,2))
plt.tight_layout()
plt.grid(color='gray', linestyle='dashed')
fig1.savefig('Check_undim_T_rhoh.png')

print("Plotting Check_undim_T_rhoh... Done!")


plt.close('all')
print("Plotting all ... Done!")