import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from uncertainties import unumpy
# import csv
# import pandas as pd

mad_to_sig = 1.4826

hdu = fits.open('/Users/ericliang/n1387/work_pub/measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits')
table = hdu[1].data

resol_ind = np.loadtxt('/Users/ericliang/n1387/work_pub/measurements/NGC1387-resolved_bool.txt')
resolve = resol_ind > 0.5
# resolve = np.array(table['resolve_spatial'] * table['resolve_spectral'], dtype=bool)

lines = open('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat').readlines()
dist = float(lines[14].strip()) # pc
alpha_co = float(lines[18].strip())

distance = table['r_gal'][resolve]

radius = table['RADRMS_EXTRAP_DECONV'][resolve]
radius_uc = table['RADRMS_EXTRAP_DECONV_uc'][resolve] # * mad_to_sig
radius_err = radius_uc * radius

mass = table['lum_extrap'][resolve] * alpha_co
mass_uc = table['mass_extrap_uc'][resolve] # * mad_to_sig
mass_err = mass_uc * mass

sigma = table['VRMS_EXTRAP_DECONV'][resolve]
sigma_uc = table['VRMS_EXTRAP_DECONV_UC'][resolve] # * mad_to_sig
sigma_err = sigma_uc * sigma

density = mass / (np.pi * radius**2)
density_err = density * np.sqrt( (mass_uc)**2 + (2*radius_uc)**2)

def pc2ang(val):
    return val/(4.84*dist/1e6)
def ang2pc(val):
    return val*4.84*dist/1e6


def radial_bin(boundary_radial, quantity, distance, error=None, weight=None, cov_factor=1.0):

    rad_low = boundary_radial[0]; rad_high = boundary_radial[-1]

    quantity_rad = np.zeros(len(boundary_radial)-1) * np.nan
    error_rad = np.zeros(len(boundary_radial)-1) * np.nan
    std_rad = np.zeros(len(boundary_radial)-1) * np.nan

    quantity_rad_all = [[] for _ in range(len(boundary_radial)-1)]
    quantityerr_rad_all = [[] for _ in range(len(boundary_radial)-1)]
    weight_all = [[] for _ in range(len(boundary_radial)-1)]

    for i,v in enumerate(quantity):
        if (distance[i] >= rad_high) or (distance[i] <= rad_low):
            continue
        if np.isnan(v):
            continue
        if (error is not None) and (np.isnan(error[i])):
            continue
        if (error is not None) and (error[i] >= v):
            continue

        index_this = np.searchsorted(boundary_radial, distance[i]) -1
        quantity_rad_all[index_this].append(v)
        if error is not None:
            quantityerr_rad_all[index_this].append(error[i])
        else:
            quantityerr_rad_all[index_this].append(np.nan)
        if weight is not None:
            weight_all[index_this].append(weight[i])
        else:
            weight_all[index_this].append(1.0)

    for i2, v2 in enumerate(quantity_rad_all):
        v2_array = np.array(v2)
        error_this = np.array(quantityerr_rad_all[i2])
        weight_this = np.array(weight_all[i2])

        number_this = len(v2_array)
        if number_this > 2:
            quantity_rad[i2] = np.nansum(weight_this * v2_array) / np.nansum(weight_this)
            error_rad[i2] = cov_factor * np.sqrt( np.nansum( weight_this**2 * error_this**2) ) / np.nansum(weight_this)
            std_rad[i2] = np.sqrt( np.nansum(weight_this * (v2_array - quantity_rad[i2] )**2) / (np.nansum(weight_this) * (number_this-1)/number_this ))
        else:
            quantity_rad[i2] = np.nan; error_rad[i2] = np.nan; std_rad[i2] = np.nan

    return quantity_rad, error_rad, std_rad    


''' plot the radial profiles '''

boundary_radial = np.linspace(0,1000,11)
plot_radial = (boundary_radial[1:]+boundary_radial[:-1]) / 2.

fig,ax= plt.subplots(2,2,figsize=(7,7))

average_this, error_this, std_this = radial_bin(boundary_radial, radius, distance, radius_err)
plt.sca(ax[0,0])
plt.scatter(pc2ang(distance), radius, s=1)
plt.errorbar(pc2ang(plot_radial), average_this, yerr=error_this, capsize=3, color='red')
plt.fill_between(pc2ang(plot_radial), average_this+std_this, average_this-std_this, alpha=0.2, color='m')
plt.ylabel(r'$R_{\rm c}$ (pc)')
plt.ylim(bottom=0)

average_this, error_this, std_this = radial_bin(boundary_radial, mass, distance, mass_err)
plt.sca(ax[1,0])
plt.scatter(pc2ang(distance), np.log10(mass), s=1)
err_asym=[np.log10(average_this/(average_this-error_this)),np.log10((average_this+error_this)/average_this)]
plt.errorbar(pc2ang(plot_radial), np.log10(average_this), yerr=err_asym, capsize=3, color='red')
plt.fill_between(pc2ang(plot_radial), np.log10(average_this+std_this), np.log10(average_this-std_this), alpha=0.2, color='m')
plt.ylabel(r'log$(M_{\rm gas}/$M$_\odot$)')

average_this, error_this, std_this = radial_bin(boundary_radial, sigma, distance, sigma_err)
plt.sca(ax[0,1])
plt.scatter(pc2ang(distance), sigma, s=1)
plt.errorbar(pc2ang(plot_radial), average_this, yerr=error_this, capsize=3, color='red')
plt.fill_between(pc2ang(plot_radial), average_this+std_this, average_this-std_this, alpha=0.2, color='m')
plt.ylabel(r'$\sigma_{\rm obs,los}$ (km s$^{-1}$)')
plt.ylim(bottom=0)

average_this, error_this, std_this = radial_bin(boundary_radial, density, distance, density_err)
plt.sca(ax[1,1])
plt.scatter(pc2ang(distance), np.log10(density), s=1)
err_asym=[np.log10(average_this/(average_this-error_this)),np.log10((average_this+error_this)/average_this)]
plt.errorbar(pc2ang(plot_radial), np.log10(average_this), yerr=err_asym, capsize=3, color='red')
plt.fill_between(pc2ang(plot_radial), np.log10(average_this+std_this), np.log10(average_this-std_this), alpha=0.2, color='m')
plt.ylabel(r'log($\Sigma_{\rm gas}$ / M$_\odot$pc$^{-2}$)')

for v1 in ax:
    for v2 in v1:
        plt.sca(v2)
        plt.xlabel(r'$R_{\mathrm{gal}} $ (arcsec)')
        plt.axvline(pc2ang(300),ls='--',color='grey')
        plt.axvline(pc2ang(600),ls='--',color='grey')

        v2.tick_params(axis='x', which='both', top=False )
        secax2 = v2.secondary_xaxis('top', functions=(ang2pc,pc2ang))
        secax2.set_xlabel(r'$R_{\rm gal}$ (pc)')

plt.subplots_adjust(hspace=0.4,wspace=0.3)

plt.savefig('/Users/ericliang/n1387/work_pub/plot/gradient.pdf')
plt.show()


''' distribution comparisons with N4526 and N4429'''
file = '/Users/ericliang/Desktop/WISDOM/literature-GMC/WISDOM_paper-and-co./2015-CARMA-N4526-tab1.txt'
data4526 = np.genfromtxt(file, skip_header=35, delimiter=[3,3, 3, 4, 2,3,5,6, 6, 6,6, 6, 6, 5,5,5,4,5,6,5 ])
resolve = ~np.isnan(data4526[:,8])
np.mean(data4526[:,8][resolve]) # size
np.std(data4526[:,8][resolve])
np.mean(data4526[:,10][resolve]) # linewidth
np.std(data4526[:,10][resolve])
np.mean(np.log10(data4526[:,14][resolve]*1e6)) # mass
np.std(np.log10(data4526[:,14][resolve]*1e6))
np.mean(np.log10( data4526[:,14][resolve]*1e6 / (np.pi * data4526[:,8][resolve]**2) )) # mass density
np.std(np.log10( data4526[:,14][resolve]*1e6 / (np.pi * data4526[:,8][resolve]**2) ))


file_4429 = '/Users/ericliang/Desktop/WISDOM/literature-GMC/WISDOM_paper-and-co./WISDOM_N4429-published/Table1.csv'
data4429= np.genfromtxt(file_4429, delimiter=',',skip_header=2)

resolve = (data4429[:,4] > 0.) & (data4429[:,6] > 0.)
print(np.mean(data4429[:,4][resolve])) # size
print(np.std(data4429[:,4][resolve]))
print(np.mean(data4429[:,6][resolve])) # linewidth
print(np.std(data4429[:,6][resolve]))
print(np.mean(np.log10(data4429[:,12][resolve]*1e5))) # mass
print(np.std(np.log10(data4429[:,12][resolve]*1e5)))
print(np.mean(np.log10( data4429[:,12][resolve]*1e5 / (np.pi * data4429[:,4][resolve]**2) ))) # mass density
print(np.std(np.log10( data4429[:,12][resolve]*1e5 / (np.pi * data4429[:,4][resolve]**2) )))


'''
# double checking
a = radius[(distance>900)&(distance<1000)]
b = radius_err[(distance>900)&(distance<1000)]

a2 = mass[(distance>900)&(distance<1000)]
b2 = mass_err[(distance>900)&(distance<1000)]

a3 = density[(distance>900)&(distance<1000)]
b3 = density_err[(distance>900)&(distance<1000)]

array_test = unumpy.uarray(a2, b2)
average_test = np.nanmean(array_test)
'''



# Global properties

# Omega
# A
# Q
# Sigma
# sigma
# M*

radial_grid, omega_grid = np.loadtxt('/Users/liangf/work/measurements_publication/profile_omega.txt',unpack=True)
galdensity_grid, galdensity_grid_err = np.loadtxt('/Users/liangf/work/measurements_publication/profile_density.txt',unpack=True,usecols=(1,2))
dispersion_grid, dispersion_grid_err = np.loadtxt('/Users/liangf/work/measurements_publication/profile_dispersion.txt',unpack=True,usecols=(1,2))
toomre_grid, toomre_grid_err = np.loadtxt('/Users/liangf/work/measurements_publication/profile_toomre.txt',unpack=True,usecols=(1,2))
oort_grid = np.loadtxt('/Users/liangf/work/measurements_publication/profile_oort.txt',usecols=(1,))
star_grid = np.loadtxt('/Users/liangf/work/measurements_publication/profile_star.txt',usecols=(1,))

boundary_radial = np.concatenate([[0], (radial_grid[1:]+radial_grid[:-1])/2., [radial_grid[-1]+np.diff(radial_grid)[-1]/2]])
mass_grid, mass_grid_err, std_this = radial_bin(boundary_radial, mass, distance, mass_err)
radius_grid, radius_grid_err, std_this = radial_bin(boundary_radial, radius, distance, radius_err)
sigma_grid, sigma_grid_err, std_this = radial_bin(boundary_radial, sigma, distance, sigma_err)
density_grid, density_grid_err, std_this = radial_bin(boundary_radial, density, distance, density_err)

cloud_prop = [ np.log10(mass_grid), radius_grid, sigma_grid, np.log10(density_grid)]
cloud_prop_err_low = [ np.log10(mass_grid/(mass_grid-mass_grid_err )), radius_grid_err, sigma_grid_err, np.log10(density_grid/(density_grid-density_grid_err))]
cloud_prop_err_up = [ np.log10((mass_grid+mass_grid_err)/mass_grid), radius_grid_err, sigma_grid_err, np.log10((density_grid+density_grid_err)/density_grid)]
gal_prop = [galdensity_grid, (omega_grid),  dispersion_grid, toomre_grid, oort_grid, np.log10(star_grid)]
gal_prop_err = [ galdensity_grid_err, None,  dispersion_grid_err, toomre_grid_err, None , None]


import matplotlib.cm as cm
from matplotlib.colors import Normalize
cmap = cm.viridis
norm = Normalize(vmin=radial_grid.min(), vmax=radial_grid.max())


# cloud_prop_label = [r'$\log(M_{\rm c}$ / M$_{\odot}$)', r'$R_{\rm c}$ (pc)',  \
#                     r'$\sigma_{\rm c}$ (km s$^{-1}$)', r'log($\Sigma_{\rm c}$ / (M$_{\odot}$ pc$^{-2}$))']
# gal_prop_label = [r'$\Sigma_{\rm mol}$ (M$_{\odot}$ pc$^{-2}$)', r'$\Omega$ (km s$^{-1}$ pc$^{-1}$)', \
#                   r'$\sigma_{\rm disc}$ (km s$^{-1}$)', r'$Q$', r'$A$ (km s$^{-1}$ pc$^{-1}$)', r'log($M_*$ / M$_{\odot}$)']


from matplotlib import rc
from matplotlib import rcParams
matplotlib.rc('text', usetex = True)
params = {'text.latex.preamble': [r'\usepackage{amsmath}']}   
plt.rcParams.update(params)

cloud_prop_label = [r'$\log(M_{\rm c}$ / M$_{\odot}$)', r'$R_{\rm c}$ (pc)',  \
                    r'$\mathit{\sigma}_{\rm c}$ (km s$^{-1}$)', r'log($\mathit{\Sigma}_{\rm c}$ / (M$_{\odot}$ pc$^{-2}$))']
gal_prop_label = [r'$\mathit{\Sigma}_{\rm mol}$ (M$_{\odot}$ pc$^{-2}$)', r'$\mathit{\Omega}$ (km s$^{-1}$ pc$^{-1}$)', \
                  r'$\mathit{\sigma}_{\rm disc}$ (km s$^{-1}$)', r'$Q$', r'$A$ (km s$^{-1}$ pc$^{-1}$)', r'log($M_*$ / M$_{\odot}$)']


fig, ax = plt.subplots(4,6,figsize=(15,8))
for i,v in enumerate(cloud_prop):
    for j,v2 in enumerate(gal_prop):
        
        plt.sca(ax[i,j])
        plt.errorbar(v2, v, xerr=gal_prop_err[j], yerr=[cloud_prop_err_low[i],cloud_prop_err_up[i]], ls=':', ecolor=cmap(norm(radial_grid)))
        pcm = plt.scatter(v2, v, s=20, c=radial_grid)

        if i == len(cloud_prop)-1:
            plt.xlabel(gal_prop_label[j])
        elif i == 0:
            ax[i,j].tick_params(labeltop=True,labelbottom=False)
        else:
            ax[i,j].axes.xaxis.set_ticklabels([])

        if j == 0:
            plt.ylabel(cloud_prop_label[i])
        elif j == len(gal_prop)-1:
            ax[i,j].tick_params(labelright=True,labelleft=False)
        else:
            ax[i,j].axes.yaxis.set_ticklabels([])

        if i == 0:
            plt.setp(ax[i,j].get_yticklabels()[1], visible=False)

# plt.subplots_adjust(wspace=0, hspace=0, right=0.95)
# cbar_ax = fig.add_axes([0.95, 0.0, 0.05, 1.0])
plt.subplots_adjust(wspace=0, hspace=0)
cb = fig.colorbar(pcm, ax=ax.ravel().tolist(), aspect=50)
cb.set_label(r'$R_{\rm gal}$ (pc)', labelpad=15, rotation=270)

plt.savefig('/Users/liangf/work/output_publication/prop_correlation-temp.pdf')
plt.show()







f=open('/Users/liangf/work/output_publication/NGC1387_gmc_table.csv','r')
reader = csv.reader(f)
labels = next(reader, None)
units = next(reader, None)
cloud_distance = []; cloud_radius = []; cloud_m = []; cloud_sigma = []; cloud_radius_err = []
for row in reader:
    if float(row[4]) > 0:
        cloud_radius.append(float(row[4]))
        cloud_radius_err.append(float(row[5]))
        cloud_sigma.append(float(row[6]))
        cloud_m.append(float(row[12]))
        cloud_distance.append(float(row[17]))
cloud_m = np.array(cloud_m) * 1e5
cloud_distance = np.array(cloud_distance); cloud_radius =np.array(cloud_radius); cloud_sigma = np.array(cloud_sigma)
cloud_Sigma = cloud_m / (np.pi * cloud_radius**2)
cloud_radius_err = np.array(cloud_radius_err)


Omega_gal_cloud = np.interp(cloud_distance, plot_radial, Omega_radial)
Sigma_gal_cloud = np.interp(cloud_distance, plot_radial, Sigma_radial)
Q_gal_cloud = np.interp(cloud_distance, plot_radial, Q)
A_gal_cloud = np.interp(cloud_distance, plot_radial, A_radial)

plt.figure()
plt.hist(cloud_distance,bins=np.sort(cloud_distance),cumulative=True);
plt.show()





cloud_prop = [np.log10(cloud_m), cloud_radius, cloud_sigma, np.log10(cloud_Sigma)]
gal_prop = [Omega_gal_cloud, Sigma_gal_cloud, Q_gal_cloud, A_gal_cloud]

fig, ax = plt.subplots(4,4,figsize=(13,7))
for i,v in enumerate(cloud_prop):
    for j,v2 in enumerate(gal_prop):
        # plt.figure()
        plt.sca(ax[i,j])
        plt.scatter(v2,v,s=2)
        plt.xlabel(gal_prop_label[j])
        plt.ylabel(cloud_prop_label[i])
        # plt.xscale('symlog')
        # plt.show()
        # break
plt.tight_layout()
plt.show()
# plt.savefig('/Users/liangf/work/prop_correlation.pdf'); plt.close()

