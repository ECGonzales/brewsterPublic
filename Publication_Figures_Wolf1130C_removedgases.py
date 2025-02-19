import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mgimg
import matplotlib.colors as colors
import scipy as sp
import numpy as np
import emcee
import testkit
import corner
import pickle as pickle
from IPython.display import display
import forwardmodel
import ciamod
import TPmod
import cloud
import band
import brewtools
from astropy.convolution import convolve, convolve_fft
from astropy.convolution import Gaussian1DKernel
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from bensconv import prism_non_uniform
from bensconv import conv_uniform_R
from bensconv import conv_uniform_FWHM
from matplotlib.lines import Line2D
import pandas as pd
from matplotlib.ticker import ScalarFormatter
import matplotlib.patheffects as pe

# ------------------------------- Load up the data and get the profile------------------
path_to_results='Results/Wolf1130C/'
figure_path='Figures/Wolf1130C/'
runname = "Wolf1130C_NC_noph3"

# error for distance and photometry for corner plot R and M
sigDist = 0.007
sigPhot = 0.02 #didn't have anything this is what Ben used for 1935

# OK finish? 1 for yes, 0 for no. #also I updated brewtools because on plieades
#I have a differen emcee version then on my laptop, or mso I though but I had to change for the next run with the .pk1 file.
fin = 1
flatendchain, flatendprobs,ndim = brewtools.get_endchain(path_to_results+runname,fin, emcee_version='3.0rc2')
print(flatendprobs)
theta_max_end = flatendchain[np.argmax(flatendprobs)]
max_end_like = np.amax(flatendprobs)
samples = flatendchain

argfile =path_to_results+runname+"_runargs.pic"

runargs = brewtools.pickle_load(argfile)
# If you're opening on a Linux box, use the code below
#with open(argfile, 'rb') as input:
#    runargs = pickle.load(input)

gases_myP,chemeq,dist, cloudtype,do_clouds,gasnum,cloudnum,inlinetemps,coarsePress,press,inwavenum,linelist,cia,ciatemps,use_disort,fwhm,obspec,proftype,do_fudge, prof,do_bff,bff_raw,ceTgrid,metscale,coscale = runargs

Tsamples = samples[:,ndim-13:]
nsamps = Tsamples.shape[0]
Tprofs = np.empty([64,Tsamples.shape[0]])
for i in range(0,nsamps):
    Tprofs[:,i] = TPmod.set_prof(1,coarsePress,press,Tsamples[i,:])

Tlays = np.empty([64,5])
for i in range(0,64):
    junk = Tprofs[i,:]
    junk2 = np.percentile(junk, [2.4,16, 50, 84,97.6],axis=0)
    junk3 = np.array(junk2)
    Tlays[i,:] = junk3[:]

# max likelihood Tprof
bestT = TPmod.set_prof(1,coarsePress,press,theta_max_end[ndim-13:])

# ---------------------------------------------------------------------------------------------------------------------
# -------------------------------- Create The Spectrum ------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
# no need to get diagnostics along with the spectrum
gnostics = 0
# Now run the model again to get your model spectrum and process to make it look like the data
shiftspec, photspec, tauspec,cfunc = testkit.modelspec(theta_max_end,runargs,gnostics)
topspec = brewtools.proc_spec(shiftspec,theta_max_end,fwhm,chemeq,gasnum,obspec)

# Now grab 500 random draws from the posterior
pltspec = np.zeros((100, obspec[0, :].size))
samp = np.empty(ndim)
sid = np.zeros(100)
for i in range(0, 100):
    sid[i] = np.random.randint(0, high=len(samples))
    samp = samples[int(sid[i]), :]
    shiftspec, photspec, tauspec, cfunc = testkit.modelspec(samp, runargs, gnostics)
    pltspec[i, :] = brewtools.proc_spec(shiftspec, samp, fwhm, chemeq, gasnum, obspec)

# get the intervals for the distribution of model spectra
specdist = np.empty([obspec[0].size, 5])
for i in range(0, obspec[0].size):
    junk = pltspec[:, i]
    junk2 = np.percentile(junk[~np.isnan(junk)], [2.4, 16, 50, 84, 97.6])
    junk3 = np.array(junk2)
    specdist[i, :] = junk3[:]


fig = plt.figure(dpi=320)
fig, axs = plt.subplots(6, gridspec_kw={'hspace': 1, 'wspace': 1})  # ,'height_ratios': [1,1,1,1,1,2]


for i in range(0, 6):
    lsize = 8

    axs[i].set(xlim=[2.8 + (i * 0.4), 3.2 + (i * 0.4)])
    # set the plot range for each plot
    xr = np.where(np.logical_and(obspec[0, :] > (2.8 + (i * 0.4)), obspec[0, :] < 3.2 + (i * 0.4)))
    xr = xr[0]

    d1, = axs[i].plot(obspec[0, xr], obspec[1, xr] / 1e-17, 'k-', label="Wolf 1130C data")
    axs[i].fill_between(obspec[0, xr], (obspec[1, xr] - obspec[2, xr]) / 1e-17, (obspec[1, xr] + obspec[2, xr]) / 1e-17,
                        facecolor='gray', alpha=0.5)
    r1, = axs[i].plot(obspec[0, xr], specdist[xr, 2] / 1e-17, 'y-', linewidth=0.8, label="Without PH3", zorder=5)
    plt.fill_between(obspec[0], specdist[:, 0], specdist[:, 4], facecolor='red', alpha=0.2)
    axs[i].fill_between(obspec[0, xr], specdist[xr, 1] / 1e-17, specdist[xr, 3] / 1e-17, facecolor='red', alpha=0.5)

    # t1, = axs[i].plot(obspec[0, xr], topspec[xr] / 1e-17, 'g-', linewidth=0.5, label="max likelihood", zorder=4)

    if i == 0:
        #         axs[i].legend(handles=[d1,r1,t1],bbox_to_anchor=(1.05,1.2),loc='upper left')
         axs[i].legend(handles=[d1, r1], bbox_to_anchor=(0.5, 1.5), loc='center', ncol=3, borderaxespad=0.2,
                       fontsize=9)


plt.ylabel(r'                                            $ F_{\lambda}$ ($10^{-17}~{\rm Wm^{-2} \mu m^{-1}}$)', fontsize=15)
plt.xlabel('Wavelength ($\mu m$)', fontsize=15)
plt.savefig(figure_path + runname +"_Pub_SPAG_SPEC_no_ph3.pdf", format='pdf', dpi=320)



# ---------------------------------------------------------------------------------------------------------------------
# -------------------------------- Create Proof of Molecules Plots ---------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
# Load up the data and models without Ph3, NH3, and H2s
data= noph3= pd.read_csv("Results/Wolf1130C/Wolf1130C_NC_noph3_observed_SPEC.dat", sep='\s+', header=None, names=["w", "f"])
winner= pd.read_csv("Results/Wolf1130C/Wolf1130C_NC_MEDIAN_SPEC.dat", sep='\s+', header=None, names=["w", "f", "err1",'err2'])
noph3= pd.read_csv("Results/Wolf1130C/Wolf1130C_NC_noph3_MEDIAN_SPEC.dat", sep='\s+', header=None, names=["w", "f", "err1",'err2'])
noh2s= pd.read_csv("Results/Wolf1130C/Wolf1130C_NC_no_h2s_MEDIAN_SPEC.dat", sep='\s+', header=None, names=["w", "f", "err1",'err2'])
nonh3= pd.read_csv("Results/Wolf1130C/Wolf1130C_NC_no_nh3_MEDIAN_SPEC.dat", sep='\s+', header=None, names=["w", "f", "err1",'err2'])


# All on one figure---- NOT so great---
fig, axs = plt.subplots(6, gridspec_kw={'hspace': 1, 'wspace': 1})  # ,'height_ratios': [1,1,1,1,1,2]

for i in range(0, 6):
    lsize = 8

    axs[i].set(xlim=[2.8 + (i * 0.4), 3.2 + (i * 0.4)])
    # set the plot range for each plot
    xr = np.where(np.logical_and(data['w']> (2.8 + (i * 0.4)), data['w'] < 3.2 + (i * 0.4)))
    xr = xr[0]

    d1, = axs[i].plot(data['w'], data['f'] / 1e-17, 'k-', linewidth=1.5,label="Wolf 1130C data")
    r1, = axs[i].plot(winner['w'], winner['f'] / 1e-17, 'y-', linewidth=1, label="Best Fit Model")
    r2, = axs[i].plot(noph3['w'], noph3['f'] / 1e-17, 'tab:green', linewidth=0.8, label="No PH3", zorder=3)
    r3, = axs[i].plot(noh2s['w'], noh2s['f'] / 1e-17, color='c', linewidth=0.5, label="No H2s", zorder=4)
    r4, = axs[i].plot(nonh3['w'], nonh3['f'] / 1e-17, color='tab:brown', linewidth=0.3, label="No NH3", zorder=5)
    # plt.fill_between(obspec[0], specdist[:, 0], specdist[:, 4], facecolor='red', alpha=0.2)
    # axs[i].fill_between(obspec[0, xr], specdist[xr, 1] / 1e-17, specdist[xr, 3] / 1e-17, facecolor='red', alpha=0.5)

    if i == 0:
        #         axs[i].legend(handles=[d1,r1,t1],bbox_to_anchor=(1.05,1.2),loc='upper left')
        axs[i].legend(handles=[d1, r1, r2, r3, r4], bbox_to_anchor=(0.5, 1.5), loc='center', ncol=5, borderaxespad=0.2,
                      fontsize=9)

axs[0].set(ylim=[0,3.4])
axs[1].set(ylim=[-0.5,3])
axs[2].set(ylim=[0,9])
axs[3].set(ylim=[0,8])
axs[4].set(ylim=[0,6.5])
axs[5].set(ylim=[0,4.8])
plt.ylabel(r'                                               $ F_{\lambda}$ ($10^{-17}~{\rm Wm^{-2} \mu m^{-1}}$)', fontsize=15)
plt.xlabel('Wavelength ($\mu m$)',fontsize=15)
#plt.savefig(runname+"_panel_spec.png",format='png', dpi=320)
plt.savefig(figure_path+"Spectra_without_gases1.png",format='png', dpi=320)


# H2S---
fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.set_size_inches(10, 6.45)
plt.gcf().subplots_adjust(bottom=0.15, left=0.15)
plt.axis([3.6, 3.76,0,6.5])
for axis in ['top', 'bottom', 'left', 'right']:  # Thicken the frame
    ax1.spines[axis].set_linewidth(1.1)

d1, = plt.plot(data['w'], data['f'] / 1e-17, 'k-', linewidth=1.5,label="Wolf 1130C data")
r1, = plt.plot(winner['w'], winner['f'] / 1e-17, 'y-', linewidth=1, label="Winning Model")
r2, = plt.plot(noh2s['w'], noh2s['f'] / 1e-17, color='c', linewidth=0.8, label="Without H2s", zorder=4)

plt.annotate('Data', xy=(3.602, 6), color='k', fontsize=20)
plt.annotate('Best Fit Model', xy=(3.602, 5.6), color='y', fontsize=20)
plt.annotate('Model without H$_2$S', xy=(3.602, 5.2), color='c', fontsize=20) #,path_effects=[pe.withStroke(linewidth=0.3, foreground="blue")]

plt.ylabel(r'$ F_{\lambda}$ ($10^{-17}~{\rm Wm^{-2} \mu m^{-1}}$)', fontsize=20)
plt.xlabel('Wavelength ($\mu m$)',fontsize=20)
ax1.tick_params(axis='both', labelsize=20, length=8, width=1.1)
plt.tight_layout()
plt.savefig(figure_path+"Spectra_without_H2S_single.png",format='png', dpi=320)

# PH3---
fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.set_size_inches(10, 6.45)
plt.gcf().subplots_adjust(bottom=0.15, left=0.15)
plt.axis([4, 4.5,2,8])
for axis in ['top', 'bottom', 'left', 'right']:  # Thicken the frame
    ax1.spines[axis].set_linewidth(1.1)

d1, = plt.plot(data['w'], data['f'] / 1e-17, 'k-', linewidth=1.5,label="Wolf 1130C data")
r1, = plt.plot(winner['w'], winner['f'] / 1e-17, 'y-', linewidth=1, label="Winning Model")
r2, = plt.plot(noph3['w'], noph3['f'] / 1e-17, color='tab:green', linewidth=0.8, label="Without Ph3", zorder=4)

plt.annotate('Data', xy=(4.34, 7.6), color='k', fontsize=20)
plt.annotate('Best Fit Model', xy=(4.34, 7.3), color='y', fontsize=20)
plt.annotate('Model without PH$_3$', xy=(4.34, 7.0), color='tab:green', fontsize=20)

plt.ylabel(r'$ F_{\lambda}$ ($10^{-17}~{\rm Wm^{-2} \mu m^{-1}}$)', fontsize=20)
plt.xlabel('Wavelength ($\mu m$)',fontsize=20)
ax1.tick_params(axis='both', labelsize=20, length=8, width=1.1)
plt.tight_layout()
plt.savefig(figure_path+"Spectra_without_Ph3_single.png",format='png', dpi=320)

# NH3---
fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.set_size_inches(10, 6.45)
plt.gcf().subplots_adjust(bottom=0.15, left=0.15)
plt.axis([2.96, 3.05,0,3])
for axis in ['top', 'bottom', 'left', 'right']:  # Thicken the frame
    ax1.spines[axis].set_linewidth(1.1)

d1, = plt.plot(data['w'], data['f'] / 1e-17, 'k-', linewidth=1.5,label="Wolf 1130C data")
r1, = plt.plot(winner['w'], winner['f'] / 1e-17, 'y-', linewidth=1, label="Winning Model")
r2, = plt.plot(nonh3['w'], nonh3['f'] / 1e-17, color='tab:brown', linewidth=0.8, label="Without NH3", zorder=4)

plt.annotate('Data', xy=(2.962, 2.8), color='k', fontsize=20)
plt.annotate('Best Fit Model', xy=(2.962, 2.6), color='y', fontsize=20)
plt.annotate('Model without NH$_3$', xy=(2.962, 2.4), color='tab:brown', fontsize=20)


plt.ylabel(r'$ F_{\lambda}$ ($10^{-17}~{\rm Wm^{-2} \mu m^{-1}}$)', fontsize=20)
plt.xlabel('Wavelength ($\mu m$)',fontsize=20)
ax1.tick_params(axis='both', labelsize=20, length=8, width=1.1)
plt.tight_layout()
plt.savefig(figure_path+"Spectra_without_NH3_single.png",format='png', dpi=320)