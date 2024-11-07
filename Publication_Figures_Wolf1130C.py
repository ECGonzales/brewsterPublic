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


#------------------------------- Load up the data and get the profile------------------
path_to_results='Results/Wolf1130C/'
figure_path='Figures/Wolf1130C/'
runname = "Wolf1130C_NC"

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
# -------------------------------- Create The PT profile ------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
#Load up some forward models for comparisons
sonoramodels = '/Users/eileengonzales/Dropbox/BDNYC/BDNYC_Research/Models/Marley2018_Sonoramodels/cmp/'
bobcatmodels = '/Users/eileengonzales/Dropbox/BDNYC/BDNYC_Research/Models/Sonora_Bobcat/SonoraBobcat_M-0.5_profiles/'
elfowlmodels = '/Users/eileengonzales/Dropbox/BDNYC/BDNYC_Research/Models/Sonora_Elf_Owl/T_dwarfs_575-1200/' #Need to download these


modP1,modT1 = np.loadtxt(sonoramodels+"t600g1000nc_m0.0.cmp",skiprows=1,usecols=(1,2),unpack=True) #solar
modP2,modT2 = np.loadtxt(bobcatmodels+"t600g1000nc_m-0.5.cmp",skiprows=1,usecols=(1,2),unpack=True) #-0.5
# modP3,modT4 = np.loadtxt(bobcatmodels+"t600g1000nc_m-0.5.cmp",skiprows=1,usecols=(1,2),unpack=True) #-0.5
modP4,modT4 = np.loadtxt("profilegrid_kz_1d7_qt_onfly_600_grav_316_mh_-1.0_cto_0.5.dat",skiprows=1,usecols=(1,2),unpack=True) #-0.5

# My usual figure configurations
fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.set_size_inches(10, 6.45)
plt.gcf().subplots_adjust(bottom=0.15, left=0.15)
# plt.axis([0., 2000.,3.0,-5.0])
plt.axis([0, 1400.,3.0,-4.0])
for axis in ['top', 'bottom', 'left', 'right']:  # Thicken the frame
    ax1.spines[axis].set_linewidth(1.1)

# ---- Plot the PT profile median and the 1 nad 2 sigma contours.--------
logP = np.log10(press)

d1, = plt.plot(Tlays[:,2],logP,'k-', label='Wolf 1130C Cloud free')
plt.fill_betweenx(logP,Tlays[:,1], Tlays[:,3], facecolor='red', alpha=0.3)
plt.fill_betweenx(logP,Tlays[:,0], Tlays[:,4], facecolor='red', alpha=0.1)

#------------ Plot the Maximum Likelihood ------
# l1, = plt.plot(bestT,logP,linestyle='dashdot', color='green',linewidth=1,label='max likelihood')

# -------- Plot the Model PT profiles ------------------------
# m1, = plt.plot(modT1,np.log10(modP1),'m-',linewidth=1,label='Sonora: 600K, logg = 5.0, [M/H]=0')
m2, = plt.plot(modT2,np.log10(modP2),'c-',linewidth=1,label='Bobcat: 600K, logg = 5.0, [M/H]=$-$0.5')
# m3, = plt.plot(modT3,np.log10(modP3),'y-',linewidth=1,label='Bobcat 500K logg = 5.5 [M/H]=+0.3')
m3, = plt.plot(modT2,np.log10(modP4),'#A04FE7',linewidth=1,label='Elf-Owl: 600K, logg = 4.5, k$_{zz}$=7, [M/H]=$-$1, C/O =0.23')

# ------- Plot the condensation curves and label -----------
mh_cc =-0.68
mns = 10.0**4/(7.45 - 0.42*logP-0.84*mh_cc)
na2s = 10.0**4/(10.05 - 0.72*logP-1.08*mh_cc)
zns = 10.0**4/(12.52 - 0.63*logP-1.26*mh_cc)
kcl = 10.0**4/(12.48 - 0.8786*logP-0.8786*mh_cc)
# nh4h2po4 = 10.0**4/(29.99-0.20*(11.0*logP + 15.0*0.0))
warm_ADP= 10**4/(20.625-0.890*logP-2.671*mh_cc) # ADP from Channon 10,000/T ≈ 20.625 − 0.890 log pt − 2.671[M/H]
h2o = 10000.0 / (40.042-3.730*logP - 3.730*mh_cc)

c1, = plt.plot(mns,logP,'--',color='#DBB93B',linewidth=1.5, label='MnS') #pick new color for this!
c2, = plt.plot(na2s,logP,'--',color='#CE96FF',linewidth=1.5,label='Na$_2$S')
c3, = plt.plot(zns,logP,'--',color='orange',linewidth=1.5, label='ZnS')
c4, = plt.plot(kcl,logP,'--',color='purple',linewidth=1.5, label='KCl')
# c7, = plt.plot(nh4h2po4,logP,'--',color='red',linewidth=1.5, label='NH$_4$H$_2$PO$_4$')
c5, = plt.plot(warm_ADP,logP,'--',color='#26940D',linewidth=1.5, label='NH$_4$H$_2$PO$_4$')
c6, = plt.plot(h2o,logP,'--',color='blue',linewidth=1.5, label='H$_2$O')

#For soloar M/H
# plt.text(1240, -1, 'MnS', rotation=-60,fontsize=14, color='#DBB93B')
# plt.text(900, -1, 'Na$_2$S', rotation=-60,fontsize=14, color='#CE96FF') ##CE96FF, #CE008F
# plt.text(750, -1, 'ZnS', rotation=-75,fontsize=14, color='orange')
# plt.text(670, -1.2, 'KCl', rotation=-65,fontsize=14, color='purple')
# plt.text(230, -0.75, 'H$_2$O', rotation=-80,fontsize=14, color='blue')

#For M/H=-0.68
plt.text(1155, -1, 'MnS', rotation=-60,fontsize=14, color='#DBB93B')
plt.text(840, -1, 'Na$_2$S', rotation=-60,fontsize=14, color='#CE96FF') ##CE96FF, #CE008F
plt.text(695, -1.2, 'ZnS', rotation=-75,fontsize=14, color='orange')
plt.text(745, 0, 'KCl', rotation=-65,fontsize=14, color='purple')
plt.text(410, -0.75, 'NH$_4$H$_2$PO$_4$', rotation=-80,fontsize=14, color='#26940D')
plt.text(215, -0.75, 'H$_2$O', rotation=-80,fontsize=14, color='blue')



#Add Legend and axes labels
plt.legend(handles=[d1,m2, m3],fontsize=14)
# plt.legend(handles=[d1,l1,m1,m2,c1,c2,c3,c4,c6],fontsize=12, loc=1)
plt.ylabel(r'log(P) (bar)', fontsize=20)
plt.xlabel('T (K)',fontsize=20)
ax1.tick_params(axis='both', labelsize=20, length=8, width=1.1)
plt.tight_layout()

plt.savefig(figure_path+runname+"_Pub_profile.pdf",format='pdf', dpi=320)


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

    d1, = axs[i].plot(obspec[0, xr], obspec[1, xr] / 1e-17, 'k-', label=runname + " data")
    axs[i].fill_between(obspec[0, xr], (obspec[1, xr] - obspec[2, xr]) / 1e-17, (obspec[1, xr] + obspec[2, xr]) / 1e-17,
                        facecolor='gray', alpha=0.5)
    r1, = axs[i].plot(obspec[0, xr], specdist[xr, 2] / 1e-17, 'y-', linewidth=0.8, label="median", zorder=5)
    plt.fill_between(obspec[0], specdist[:, 0], specdist[:, 4], facecolor='red', alpha=0.2)
    axs[i].fill_between(obspec[0, xr], specdist[xr, 1] / 1e-17, specdist[xr, 3] / 1e-17, facecolor='red', alpha=0.5)

    t1, = axs[i].plot(obspec[0, xr], topspec[xr] / 1e-17, 'g-', linewidth=0.5, label="max likelihood", zorder=4)

    if i == 0:
        #         axs[i].legend(handles=[d1,r1,t1],bbox_to_anchor=(1.05,1.2),loc='upper left')
        axs[i].legend(handles=[d1, r1, t1], bbox_to_anchor=(0.5, 1.5), loc='center', ncol=3, borderaxespad=0.2,
                      fontsize=9)


plt.ylabel(r'                                          $ F_{\lambda}$ ($10^{-17}~{\rm Wm^{-2} \mu m^{-1}}$)', fontsize=15)
plt.xlabel('Wavelength ($\mu m$)', fontsize=15)
plt.savefig(figure_path + runname +"_Pub_SPAG_SPEC.pdf", format='pdf', dpi=320)



# ---------------------------------------------------------------------------------------------------------------------
# -------------------------------- Create The VMR ------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
# This saves out the gas mixing ratios
gasVMR = np.ones([gasnum.size, 3])
for i in range(0, gasnum.size):
    gasVMR[i, :] = np.percentile(samples[:, i], [16, 50, 84])

# np.savetxt(runname+'_VMRs.dat', np.c_[gasVMR[:,0],gasVMR[:,1],gasVMR[:,2]])
np.savetxt(path_to_results + runname + '_VMRs.dat', np.c_[gasVMR[:, 0], gasVMR[:, 1], gasVMR[:, 2]])
# grab BFF and Chemical grids
gaslist = ['h2o','ch4','co','co2','nh3','h2s','ph3']
print(gaslist)
chemeq=1
bff_raw,ceTgrid,metscale,coscale,gases_myP = testkit.sort_bff_and_CE(chemeq,"chem_eq_tables_P3K.pic",press,gaslist)

#This grid looks crazy when plotting at low t ask CHannon about it.
#bff_raw,ceTgrid,metscale,coscale,gases_myP = testkit.sort_bff_and_CE(chemeq,"chem_eq_tables_grid_500_6000_P3K.pic",press,gaslist)

# what metallicity and C/O ratio do we want to plot?
mh  = -0.7
co = (0.26/0.45)

logP, profT = np.loadtxt(path_to_results+runname+"_profile.dat",unpack=True)

nlayers = press.size
mfit = interp1d(metscale,gases_myP,axis=0)
gases_myM = mfit(mh)
cfit = interp1d(coscale,gases_myM,axis=0)
invmr = cfit(co)
ng = invmr.shape[2]
ngas = len(gaslist)
ab = np.zeros([nlayers,ngas],dtype='d')
temp= profT #Tlays[:,2]
for p in range(0,nlayers):
    for g in range(3,ng):
        tfit = InterpolatedUnivariateSpline(ceTgrid,invmr[:,p,g])
        ab[p,g-3]= tfit(temp[p])

vmr = gasVMR

fig6=plt.figure(dpi=320)
# ax = fig.add_subplot(1,1,1)
plt.axis([-12, -2, 3.0, -4.0])

colors=['#00BFBF','#87405C','#FF5186', 'orange', '#EE7F9C', 'blue','#678D56']
# gas_names_list =['h2o', 'ch4', 'co', 'co2', 'nh3', 'h2s', 'ph3']
gas_names_list =['H$_{2}$O', 'CH$_{4}$', 'CO', 'CO$_{2}$',
                 'NH$_{3}$', 'H$_{2}$S', 'PH$_{3}$']


for i, (color, gas_names_list) in enumerate (zip(colors, gas_names_list)):
    plt.plot(ab[:,i],logP,'--',color=color,linewidth=1.5)
    plt.plot([vmr[i,1],vmr[i,1]],[2.4,-4.0],'-',color=color,label=gas_names_list.upper())
    plt.fill_betweenx([2.4,-4.0],[vmr[i,0],vmr[i,0]],[vmr[i,2],vmr[i,2]], facecolor=color, alpha=0.3,linewidth=0)

plt.fill_between([-12.0,12.0],np.log10(70),np.log10(0.1),facecolor='grey',alpha=0.1) # block for the photosphere
plt.annotate('photosphere', xy=(-3.6, 1.8), color='k', fontsize=8)
# plt.annotate('- - [M/H]=-0.5', xy=(-3.45, -0.2), color='k', fontsize=10)
# plt.annotate('Model', xy=(-3.3, 0.2), color='k', fontsize=10)

# Get automatic labels
handles, labels =plt.gca().get_legend_handles_labels()

#Create my Own label and append to the legend
# line = Line2D([0], [0], linestyle='--', label='[M/H]=-0.7 \n Model', color='k')
# handles.extend([line])
# plt.legend(handles=handles, loc='upper right',fontsize=8)


# Add legend and labels that are shared over all the subplots
# create lines for custom legend
custom_lines = [Line2D([0], [0], color='k', ls='-'),
                Line2D([0], [0], color='k', ls='dashed')]
legend=plt.legend(custom_lines, ['Retrieved Abundance', '[M/H]=-0.7, C/0=0.26 Model Abundance'], loc='lower left',
                  bbox_to_anchor= (0, 1.01), ncol=2, borderaxespad=0, fontsize=9)


plt.annotate('H$_\mathrm{2}$O', xy=(-3, -3.5), color='#00BFBF', fontsize=15)
plt.annotate('CH$_\mathrm{4}$', xy=(-3, -3), color='#87405C', fontsize=15)
plt.annotate('CO', xy=(-3, -2.5), color='#FF5186', fontsize=15)
plt.annotate('CO$_\mathrm{2}$', xy=(-3, -2), color='orange', fontsize=15)
plt.annotate('NH$_\mathrm{3}$', xy=(-3, -1.5), color='#EE7F9C', fontsize=15)
plt.annotate('H$_\mathrm{2}$S', xy=(-3, -1), color='blue', fontsize=15)
plt.annotate('PH$_\mathrm{3}$', xy=(-3, -0.5), color='#678D56', fontsize=15)


#Add axes Labels
plt.ylabel('log$\,$P  (bars)',fontsize=15)
plt.xlabel('log$\,f_{gas}$',fontsize=15)
plt.savefig(figure_path + runname + '_Pub_abundances.pdf', format='pdf', dpi=320, bbox_inches='tight')

#Post process corner plot
samples = brewtools.pickle_load(path_to_results+runname_post+'_postprod.pk1')
# You'll need to edit this cell to make it work for the gases you've used
# e.g. change the cut you take out of the samples array to get gases and gravity
# gassamples = samples[:,0:14].copy()
gassamples = samples[:, 0:15].copy()
# set up boundaries for histograms
rh2o = (-3.9, -3.55)
rch4 = (-4.4, -4.1)
rco = (-7.2, -6.8)
rco2 = (-12, -8.5)
rnh3 = (-6, -5)
rh2s = (-6, -4)
rph3 = (-7.2, -6.8)
rmet = (-1, -0.5)
rmet_full = (-1, -0.5)
rco_rat = (0.1, 0.4)
rlogg = (4.8, 5.5)
rRad = (0.8, 1)
rMass = (15, 100)
rvrad = (-35, -30)
rvsini = (0, 80)

# get M/H and C/O
h2o = gassamples[:, 0]
ch4 = gassamples[:, 1]
co = gassamples[:, 2]
co2 = gassamples[:, 3]
nh3 = gassamples[:, 4]
h2s = gassamples[:, 5]
ph3 = gassamples[:, 6]

# first get the C/O

O = 10 ** h2o + 10 ** co + 2. * 10 ** (co2)
C = 10 ** (co) + 10 ** (co2) + 10 ** (ch4)

CO_ratio = C / O

# won't bother with rest of the elements as they're not fully accounted for
# i.e. N is hiding in N2, S maybe hiding elsewhere,
# and we haven't detected PH3

# Eileen will here for Wolf1130C to compare to value.
N = 10 ** nh3
S = 10 ** h2s
P = 10 ** ph3

# Determine "fraction" of H2 in the L dwarf
gas_sum = 10 ** h2o + 10 ** co + 10 ** co2 + 10 ** ch4  # +10**ph3 + 10**h2s + 10**nh3
fH2 = (1 - gas_sum) * 0.84  # fH2/(fH2+FHe) = 0.84
fH = 2. * fH2

# Determine linear solar abundance sum of elements in our T dwarf
# abundances taken from Asplund+ 2009
solar_H = 12.00
solar_O = 10 ** (8.69 - solar_H)
solar_C = 10 ** (8.43 - solar_H)
# solar_Ti = 10**(4.95-solar_H)
# solar_V = 10**(3.93-solar_H)
# solar_Cr = 10**(5.64-solar_H)
# solar_Fe = 10**(7.50-solar_H)
# solar_NaK = 10**(6.24-solar_H) + 10**(5.03-solar_H)
solar_N = 10 ** (7.83 - solar_H)
solar_S = 10 ** (7.12 - solar_H)
solar_P = 10 ** (5.41 - solar_H)

# Calculate the metallicity fraction in the star and the same for the sun and then make the ratio
# Eileen will include Ph3 as it is constrained!
metallicity_target = (O / fH) + (C / fH) + (P / fH)
metallicity_sun = solar_O + solar_C + solar_P

MH = np.log10(metallicity_target / metallicity_sun)

# Calculate the metallicity fraction in the star and the same for the sun using all retrieved metals and then make the ratio
metallicity_target_full = (O / fH) + (C / fH) + (P / fH) + (N / fH) + (S / fH)
metallicity_sun_full = solar_O + solar_C + solar_P + solar_N + solar_S

MH_full = np.log10(metallicity_target_full / metallicity_sun_full)

# Get the radius and mass
r2d2 = samples[:, 8].copy()
logg = samples[:, 7].copy()
# Radius
sigR2D2 = sigPhot * r2d2 * (-1. / 2.5) * np.log(10.)

sigD = sigDist * 3.086e16
D = dist * 3.086e16

R = np.sqrt(((np.random.randn(len(samples[:, 0])) * sigR2D2) + r2d2)) * (
            (np.random.randn(len(samples[:, 0])) * sigD) + D)

# and mass
g = (10. ** logg) / 100.
M = (R ** 2 * g / (6.67E-11)) / 1.898E27

R = R / 71492e3

gassamples[:, 7] = MH
gassamples[:, 8] = CO_ratio
gassamples[:, 9] = logg
gassamples[:, 10] = R
gassamples[:, 11] = M
gassamples[:, 12:14] = samples[:, 9:11].copy()
gassamples[:, 14] = MH_full
# now make an array of the bounds to give to corner plot
bnds = [rh2o, rch4, rco, rco2, rnh3, rh2s, rph3, rmet, rco_rat, rlogg, rRad, rMass, rvrad, rvsini, rmet_full]
# range = bnds,

fig = corner.corner(gassamples, scale_hist=False, plot_datapoints=False, range=bnds,
                    labels=["H$_2$O", "CH$_4$", "CO", "CO$_2$", "NH$_3$", "H$_2$S", "PH$_3$", "[M/H]", "C/O",
                            "$\log g$",
                            "$R / R_{Jup}$", "$M / M_{Jup}$", "$v_{rad}$", "$v \sin i$", "[M/H]_full"],
                    quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_kwargs={"fontsize": 20, "loc": "left"}, label_kwargs={"fontsize": 20})
# plt.savefig(runname+"_MRgascorner.png",bbox_inches='tight',format='png', dpi=320)
plt.savefig(figure_path + runname + "_ext" + ext_num + "_gascorner_full.png", format='png', dpi=320)