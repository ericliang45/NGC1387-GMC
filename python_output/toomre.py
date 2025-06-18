import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from uncertainties import unumpy

# from scipy.stats import median_abs_deviation
# from statsmodels.stats.weightstats import DescrStatsW
# import weightedstats as ws
# import lfh

solar_mass = 1.98847e30 # kg
pc_m = 3.086e+16
g_const2 = 4.302e-3 # (pc / Msun) * (km/s)^2
# g_const = 6.67408e-11 # m^3 / kg / s^2


gal_property = np.genfromtxt('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat', comments=';',dtype=str)
pa = float(gal_property[0]) # degree
inclination = float(gal_property[1]) # degree
dist = float(gal_property[4]) # pc
co_factor = float(gal_property[6]) # M_solar / pc^2 / (K km/s)

moments_root = '/Users/ericliang/n1387/work_pub/plot/moment_maps/rms_only-202406/'

hdr = fits.getheader(moments_root+'NGC1387_mom0.fits')
Intensity_map = fits.getdata(moments_root+'/NGC1387_mom0.fits')
Intensity_error = fits.getdata(moments_root+'/NGC1387_mom0-err9.fits')

# Temporary, testing new masks
# Intensity_map_old = fits.getdata('/Users/ericliang/n1387/work_pub/plot/moment_maps/fits//NGC1387_mom0.fits')
# co_factor = float(gal_property[6]) * np.nansum(Intensity_map_old) / np.nansum(Intensity_map)
# Temporary above

sigma_map = fits.getdata(moments_root+'/NGC1387_mom2.fits')
sigma_error = fits.getdata(moments_root+'/NGC1387_mom2-err8.fits')
shape = sigma_map.shape
center = [int(len(sigma_map)/2), int(len(sigma_map[0])/2)]
bmaj_value = hdr['BMAJ'] * 3600 # arcsec
bmin_value = hdr['BMIN'] * 3600 # arcsec
cellsize = hdr['CDELT2'] * 3600 # arcsec per spaxel

# plt.figure()
# plt.imshow(Sigma_map/Sigma_error_ori)
# plt.colorbar()
# plt.title('Significance')
# plt.show()

''' for other beam conversion'''
# from astropy import units as u
# freq =  230.538 *u.GHz
# bmaj = 7.2 * u.arcsec; bmin = 4.6 * u.arcsec
# beam_area = 1.1330900354567983 * (bmaj*bmin)
# factor = (1.0 * (u.Jy/beam_area).to(u.K, equivalencies=u.brightness_temperature(freq)) ).value 
# print(factor)
''' end '''

rest_freq = 230.538 # GHz

freq = rest_freq*u.GHz
bmaj = bmaj_value * u.arcsec; bmin = bmin_value * u.arcsec
beam_area = 1.1330900354567983 * (bmaj*bmin)
factor = (1.0 * (u.Jy/beam_area).to(u.K, equivalencies=u.brightness_temperature(freq)) ).value 
# print(0.41 * factor) # mJy * factor = mK

tb_map = Intensity_map * factor
Sigma_map = tb_map * co_factor * np.cos(np.deg2rad(inclination)) # solar mass per pc^2
Sigma_error = Intensity_error * factor * co_factor * np.cos(np.deg2rad(inclination))



x = np.arange(0,shape[0],1); y = np.arange(0,shape[1],1)
xv, yv = np.meshgrid(x, y)
delX = xv - center[1]
delY = yv - center[0]
ang2Rot = 90 - (pa -180)     # in [degree]
c = np.cos(np.radians(ang2Rot))
s = np.sin(np.radians(ang2Rot))
majComp = delX * c - delY * s
minComp = delX * s + delY * c
minProj = minComp/np.cos(np.radians(0))
r_flat_spaxel = np.sqrt((majComp * majComp) + (minProj * minProj))
r_flat_arcsec = r_flat_spaxel * cellsize
pc_per_arcsec = 1. / 3600. / 180. * np.pi * dist # pc
r_flat_pc_noi = r_flat_arcsec * pc_per_arcsec # no inclination correction, for stellar mass profile production


x = np.arange(0,shape[0],1); y = np.arange(0,shape[1],1)
xv, yv = np.meshgrid(x, y)
delX = xv - center[1]
delY = yv - center[0]
ang2Rot = 90 - (pa -180)     # in [degree]
c = np.cos(np.radians(ang2Rot))
s = np.sin(np.radians(ang2Rot))
majComp = delX * c - delY * s
minComp = delX * s + delY * c
minProj = minComp/np.cos(np.radians(inclination))
r_flat_spaxel = np.sqrt((majComp * majComp) + (minProj * minProj))
r_flat_arcsec = r_flat_spaxel * cellsize
pc_per_arcsec = 1. / 3600. / 180. * np.pi * dist # pc
r_flat_pc = r_flat_arcsec * pc_per_arcsec # 



# plt.figure()
# plt.title('Galactocentri distance')
# plt.imshow(r_flat_arcsec);
# cb = plt.colorbar()
# cb.set_label('arcsec')
# plt.show()

# plt.figure()
# plt.title('Galactocentri distance')
# plt.imshow(r_flat_pc)
# cb=plt.colorbar()
# cb.set_label('pc')
# plt.show()


''' Checking center '''
# kernel = Gaussian2DKernel(x_stddev=2)
# smooth_cloud = convolve(Sigma_map, kernel)
# plt.figure()
# plt.imshow(r_flat_pc,vmin=0,vmax=10)
# cb=plt.colorbar()
# cb.set_label('pc')
# plt.contour(smooth_cloud,levels=[7.5e-3,3.0e-2, 6.5e-2, 0.1],cmap='Greys',alpha=1.0) #  transform=ax.get_transform(WCS(cloud_hdu[0].header)
# plt.show()



def radial(map_distance, boundary_radial, map_quantity, map_error=None, map_weight=None, cov_factor=1., value_min=0.):
    rad_low = boundary_radial[0]; rad_high = boundary_radial[-1]

    quantity_rad = np.zeros(len(boundary_radial)-1) * np.nan
    error_rad = np.zeros(len(boundary_radial)-1) * np.nan
    std_rad = np.zeros(len(boundary_radial)-1) * np.nan

    quantity_rad_all = [[] for _ in range(len(boundary_radial)-1)]
    quantityerr_rad_all = [[] for _ in range(len(boundary_radial)-1)]
    weight_all = [[] for _ in range(len(boundary_radial)-1)]

    for i,v1 in enumerate(map_quantity):
        for j, v in enumerate(v1):
            if (map_distance[i,j] >= rad_high) or (map_distance[i,j] <= rad_low):
                continue
            if v <= value_min:
                continue
            if np.isnan(v):
                continue
            if (map_error is not None) and (np.isnan(map_error[i,j])):
                continue

            index_this = np.searchsorted(boundary_radial, map_distance[i,j]) -1
            quantity_rad_all[index_this].append(v)
            if map_error is not None:
                quantityerr_rad_all[index_this].append(map_error[i,j])
            else:
                quantityerr_rad_all[index_this].append(np.nan)
            if map_weight is not None:
                weight_all[index_this].append(map_weight[i,j])
            else:
                weight_all[index_this].append(1.0)

    for i2, v2 in enumerate(quantity_rad_all):
        v2_array = np.array(v2)
        error_this = np.array(quantityerr_rad_all[i2])
        weight_this = np.array(weight_all[i2])

        number_this = len(v2_array)
        quantity_rad[i2] = np.nansum(weight_this * v2_array) / np.nansum(weight_this) # weighted mean
        error_rad[i2] = cov_factor * np.sqrt( np.nansum( weight_this**2 * error_this**2) ) / np.nansum(weight_this)

        # Alternative, exactly the same
        # array_err = unumpy.uarray( v2_array, error_this )
        # result = np.sum(weight_this * array_err) / np.sum(weight_this)
        # quantity_rad[i2] = unumpy.nominal_values(result)
        # error_rad[i2] = unumpy.std_devs(result) * cov_factor

        if number_this > 1:
            std_rad[i2] = np.sqrt( np.nansum(weight_this * (v2_array - quantity_rad[i2] )**2) / (np.nansum(weight_this) * (number_this-1)/number_this ))

    return quantity_rad, error_rad, std_rad


def ang2pctrans(val): # angle in arcsec, dist in pc
    return val*4.84*dist/1e6

def ang2pctrans_inv(val): # angle in arcsec, dist in pc
    return val/(4.84*dist/1e6)


boundary_radial = np.linspace(0,1150,116)
plot_radial = (boundary_radial[1:]+boundary_radial[:-1]) / 2.


''' check there's no nan; nan should be zero '''
# plt.figure()
# plt.imshow(Sigma_map)
# plt.show()

# plt.figure()
# plt.imshow(sigma_map)
# plt.show()



''' Sigma radial profile '''

Sigma_radial, Sigma_radial_error, Sigma_radial_std = radial(r_flat_pc, boundary_radial, Sigma_map, Sigma_error, cov_factor=np.sqrt(beam_area.value / cellsize**2)) # mass surface density # , map_weight= 1./Sigma_error**2

# check the distribution is ordinary
# plt.figure()
# plt.scatter(r_flat_pc, Sigma_map,s=1)
# plt.show()

plt.figure(figsize=(4,3))
plt.fill_between(ang2pctrans_inv(plot_radial), Sigma_radial+Sigma_radial_std, Sigma_radial-Sigma_radial_std, alpha=0.3)
plt.errorbar(ang2pctrans_inv(plot_radial), Sigma_radial, yerr=Sigma_radial_error, capsize=3,) # ecolor='black'
# plt.scatter(plot_radial, Sigma_radial, s=2)
plt.axvline(ang2pctrans_inv(300),ls='--',color='gray');plt.axvline(ang2pctrans_inv(600),ls='--',color='gray')
# plt.axvline(525,ls='--',color='gray');plt.axvline(600,ls='--',color='gray')
# plt.axvline(700,ls='--',color='gray');plt.axvline(600,ls='--',color='gray')
plt.xlabel(r'$R_{\mathrm{gal}} $ (arcsec)');plt.ylabel(r'$\Sigma_{\mathrm{mol}}$ (M$_\odot$ pc$^{-2}$)')
plt.xlim(left=-0.2,right=11);plt.ylim(bottom=0) # some emission outside 11", but not reliable
ax1 = plt.gca()
ax1.tick_params(axis='x', which='both', top=False )
secax2 = ax1.secondary_xaxis('top', functions=(ang2pctrans,ang2pctrans_inv))
secax2.set_xlabel(r'$R_{\mathrm{gal}} $ (pc)')
# plt.savefig('/Users/ericliang/n1387/work_pub/plot/density_profile.pdf')
plt.show()



boundary_radial = np.linspace(0,1000,81)
plot_radial = (boundary_radial[1:]+boundary_radial[:-1]) / 2.
Sigma_radial, Sigma_radial_error, Sigma_radial_std = radial(r_flat_pc, boundary_radial, Sigma_map, Sigma_error, cov_factor=np.sqrt(beam_area.value / cellsize**2)) # mass surface density # , map_weight= 1./Sigma_error**2
sigma_radial, sigma_radial_error, sigma_radial_std = radial(r_flat_pc, boundary_radial, sigma_map, sigma_error, map_weight=Sigma_map, cov_factor=np.sqrt(beam_area.value / cellsize**2)) # velocity dispersion     


plt.figure(figsize=(4,3))
plt.fill_between(ang2pctrans_inv(plot_radial), sigma_radial+sigma_radial_std, sigma_radial-sigma_radial_std, alpha=0.3)
plt.errorbar(ang2pctrans_inv(plot_radial), sigma_radial, yerr=sigma_radial_error, capsize=3,) # ecolor='black'
plt.axvline(ang2pctrans_inv(300),ls='--',color='gray');plt.axvline(ang2pctrans_inv(600),ls='--',color='gray')
plt.xlabel(r'$R_{\mathrm{gal}} $ (arcsec)');plt.ylabel(r'$\sigma$ (km s$^{-1}$)')
ax1 = plt.gca()
ax1.tick_params(axis='x', which='both', top=False )
secax2 = ax1.secondary_xaxis('top', functions=(ang2pctrans, ang2pctrans_inv))
secax2.set_xlabel(r'$R_{\mathrm{gal}} $ (pc)')
# plt.savefig('/Users/ericliang/n1387/work_pub/plot/vdis.pdf')
plt.show()



curve_rad, curve_v = np.genfromtxt('/Users/ericliang/n1387/work_pub/data/NGC1387_velocity_curve.dat',comments=';',delimiter=',',unpack=True)

Omega = curve_v / curve_rad # km / s / pc
curve_rad_clip = (curve_rad[1:]+curve_rad[:-1]) / 2.

# new trial

gal_grid = np.linspace(1,1500,800)
curve_v_inter = np.interp(gal_grid, curve_rad, curve_v)
gal_grid_clip = (gal_grid[1:]+gal_grid[:-1]) / 2.
curve_v_inter_clip = np.interp(gal_grid_clip,gal_grid,curve_v_inter)
Omega_inter = curve_v_inter / gal_grid # km / s / pc
Omega_inter_clip = (Omega_inter[:-1] + Omega_inter[1:]) / 2.
dO_dR_inter = np.diff(Omega_inter) / np.diff(gal_grid)
dO2_dR_inter = np.diff(Omega_inter**2) / np.diff(gal_grid)

A_inter =  - gal_grid_clip / 2. * dO_dR_inter # 
A_radial = np.interp(plot_radial,gal_grid_clip,A_inter)

kappa = np.sqrt( gal_grid_clip * dO2_dR_inter + 4. * Omega_inter_clip ** 2 )
kappa_radial = np.interp(plot_radial, gal_grid_clip, kappa)

# old codes
# Omega_clip = (Omega[:-1] + Omega[1:]) / 2.
# dO2_dR = np.diff(Omega**2) / np.diff(curve_rad)

# kappa = np.sqrt( curve_rad_clip * dO2_dR + 4. * Omega_clip ** 2 )

# Omega_radial = np.interp(plot_radial, curve_rad_clip, Omega_clip)
# kappa_radial = np.interp(plot_radial, curve_rad_clip, kappa)

# another formula for kappa, Pringle J. E., King A., 2007, Astrophysical Flows p.161
beta = np.diff( np.log(curve_v_inter) ) / np.diff( np.log(gal_grid) )
kappa_2 = np.sqrt(2) * curve_v_inter_clip / gal_grid_clip * np.sqrt( 1+beta )
kappa_2_radial = np.interp(plot_radial, gal_grid_clip, kappa_2)

plt.figure()
plt.scatter(kappa_radial, kappa_2_radial)
plt.xlabel('kappa classic')
plt.ylabel('kappa Pringle and King 2007')
plt.plot([0.4,68],[0.4,68],ls='--')
plt.show()


plt.figure()
plt.plot(curve_rad / dist /np.pi*180.*3600,curve_v)
plt.gca().set_ylim([0.,400.])
plt.gca().set_xlim([0.,10.])
plt.xlabel('R [arcsec]')
plt.ylabel('v [km/s]')
plt.show()

plt.figure()
plt.plot(curve_rad / dist /np.pi*180.*3600,Omega)
plt.gca().set_ylim([0.1,400.])
plt.gca().set_xlim([0.,10.])
plt.xlabel('R [arcsec]')
plt.ylabel(r'$\Omega$ [km/s/pc]')
plt.yscale('log')
plt.show()


plt.figure()
plt.plot(plot_radial,A_radial)
plt.gca().set_ylim([0.1,300.])
plt.gca().set_xlim([0.,1000.])
plt.xlabel('R (pc)')
plt.ylabel(r'Oort $A$ (km/s/pc)')
plt.yscale('log')
plt.show()

''' extra, timescale estimate '''
# # eq 9 of N404 paper
# km_per_pc = 3.086e+13
# shear_time_myr = 1./(2*A_radial/km_per_pc) / (1e6 * 31536000) # 1.5-2.5 Myr

# # g_const3 = g_const2 * ...
# # above eq 14 of N404 paper & "1-12 dynamical timescale" A&A 510, A110 (2010) # dyn_time = np.sqrt(3*np.pi/32/g_const3/rho)
# dyn_time2 = 2*20*km_per_pc / 4 / (1e6 * 31536000) # 5-10 Myr # 

# # eq 6 of N404 paper
# collision_time = A_radial**(1/3) * 10**5.5**(1/3) / (2**(1/3)*g_const2**(2/3)*Sigma_radial) * km_per_pc / (1e6 * 31536000) # 5-10 Myr
''' extra, timescale estimate end '''

# plt.figure()
# plt.plot(gal_grid_clip, kappa) # curve_rad_clip
# x_low = 100; x_high = 1000
# plt.xlim(x_low,x_high)
# index_low = np.min(np.where(gal_grid_clip>x_low)[0]); index_high = np.max(np.where(gal_grid_clip<x_high)[0])
# ylow = np.min(kappa[index_low:index_high]); yhigh=np.max( kappa[index_low:index_high] ); ydiff=yhigh - ylow
# plt.ylim(ylow-0.2*ydiff, yhigh+0.2*yhigh)
# plt.show()

plt.figure()
plt.plot(plot_radial, kappa_radial)
plt.ylim(0,5)
plt.xlim(0,1050)
plt.xlabel('Galactocentric R [pc]')
plt.ylabel(r'$\kappa$ [km/s/pc]')
# plt.savefig('/Users/ericliang/n1387/work_pub/plot/kappa.pdf')
plt.show()


# Q = kappa_radial * sigma_radial / np.pi / g_const2 / Sigma_radial
# Q2 = kappa_radial * sigma_radial / np.pi / g_const / ( Sigma_radial * solar_mass) * 1e3 * 1e3 * pc_m  # tested, equivalent to Q

sigma_radial_u = unumpy.uarray(sigma_radial, sigma_radial_error)
Sigma_radial_u = unumpy.uarray(Sigma_radial, Sigma_radial_error)
Qu = kappa_radial * sigma_radial_u / np.pi / g_const2 / Sigma_radial_u
Q = unumpy.nominal_values(Qu)
Q_err = unumpy.std_devs(Qu)

Q_err_new = np.sqrt( Q_err**2 + (Q * 0.3)**2 )


# Sigma_radial[Sigma_radial==0.] =np.nan

plt.figure()
plt.plot(plot_radial, Q)
# plt.errorbar(plot_radial, Q, yerr=unumpy.std_devs(Q_un))
plt.xlabel(r'R$_{\rm gal}$ (pc)')
plt.ylabel(r'$Q$')
plt.ylim(0,4)
plt.axhline(1,ls='--',color='gray')
plt.show()


plt.figure(figsize=(4*1.2,3*1.2))
plt.fill_between(ang2pctrans_inv(plot_radial[1:]), Q[1:]+Q_err_new[1:], Q[1:]-Q_err_new[1:], alpha=0.3, color='gray')
plt.errorbar(ang2pctrans_inv(plot_radial[1:]), Q[1:], yerr=Q_err[1:], capsize=3, label=r'Toomre parameter $Q$') # ecolor='black'
plt.axvline(ang2pctrans_inv(300),ls='--',color='k');plt.axvline(ang2pctrans_inv(600),ls='--',color='k') # ,label='Region boundary'
plt.xlabel(r'$R_{\mathrm{gal}} $ (arcsec)');plt.ylabel(r'$Q$')
ax1 = plt.gca()
ax1.tick_params(axis='x', which='both', top=False )
secax2 = ax1.secondary_xaxis('top', functions=(ang2pctrans, ang2pctrans_inv))
secax2.set_xlabel(r'$R_{\mathrm{gal}} $ (pc)')
plt.ylim(0,8); plt.xlim(0,14)
plt.text(ang2pctrans_inv(1100),0.5,'Unstable',color='gray'); plt.text(ang2pctrans_inv(1100),1.2,'Stable',color='gray')
plt.axhline(1,ls='-',color='gray') # ,label='Instability threshold'
# plt.legend()
# plt.savefig('/Users/ericliang/n1387/work_pub/plot/toomre.pdf')
plt.show()




# To output Omega, A, Sigma, sigma

'''
boundary_radial = np.linspace(0,1000,11)
plot_radial = (boundary_radial[1:]+boundary_radial[:-1]) / 2.

Sigma_radial, Sigma_radial_error, Sigma_radial_std = radial(r_flat_pc, boundary_radial, Sigma_map, Sigma_error, cov_factor=np.sqrt(beam_area.value / cellsize**2)) # mass surface density # , map_weight= 1./Sigma_error**2
sigma_radial, sigma_radial_error, sigma_radial_std = radial(r_flat_pc, boundary_radial, sigma_map, sigma_error, map_weight=Sigma_map, cov_factor=np.sqrt(beam_area.value / cellsize**2)) # velocity dispersion     


Omega_output = np.interp(plot_radial, gal_grid, Omega_inter)

A_output = np.interp(plot_radial, gal_grid_clip, A_inter)

kappa_output =  np.interp(plot_radial, gal_grid_clip, kappa)

sigma_radial_u = unumpy.uarray(sigma_radial, sigma_radial_error)
Sigma_radial_u = unumpy.uarray(Sigma_radial, Sigma_radial_error)
Qu_output = kappa_output * sigma_radial_u / np.pi / g_const2 / Sigma_radial_u
Q_output = unumpy.nominal_values(Qu_output)
Q_output_err = unumpy.std_devs(Qu_output)

header = 'r_gal (pc); angular velocity (km/s/pc)'
data_this = np.array([plot_radial, Omega_output]).T
np.savetxt('/Users/liangf/work/measurements_publication/profile_omega.txt', data_this, header=header)

header = 'r_gal (pc); surface mass density (M_solar / pc^2); measurment error; standard deviation'
data_this = np.array([plot_radial, Sigma_radial, Sigma_radial_error, Sigma_radial_std]).T
np.savetxt('/Users/liangf/work/measurements_publication/profile_density.txt', data_this, header=header)

header = 'r_gal (pc); average velocity dispersion (km/s); measurment error; standard deviation'
data_this = np.array([plot_radial, sigma_radial, sigma_radial_error, sigma_radial_std]).T
np.savetxt('/Users/liangf/work/measurements_publication/profile_dispersion.txt', data_this, header=header)

header = 'r_gal (pc); Oort constant A (km/s/pc)'
data_this = np.array([plot_radial, A_output]).T
np.savetxt('/Users/liangf/work/measurements_publication/profile_oort.txt', data_this, header=header)

header = 'r_gal (pc); Toomre parameter Q; measurment error'
data_this = np.array([plot_radial, Q_output, Q_output_err]).T
np.savetxt('/Users/liangf/work/measurements_publication/profile_toomre.txt', data_this, header=header)
'''


''' stellar mass '''
mge_file = '/Users/ericliang/n1387/work_pub/data/NGC1387_mge_gauss.csv'
gal_property = np.genfromtxt('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat', comments=';',dtype=str)

# mge_file = '/Users/ericliang/n1387/work/data/NGC1387_mge_gauss-linear.csv' # for comparison
# gal_property = np.genfromtxt('/Users/ericliang/n1387/work/data/NGC1387_galpar.dat', comments=';',dtype=str)  # for comparison

pa = float(gal_property[0])
ml_par = float(gal_property[7])
ellipse = np.loadtxt(mge_file, delimiter=",", comments='#')
arcsec_per_pixel = 0.04
# dist

amap = np.zeros([500,500]) # for testing

# center = len(Sigma_map)/2
x = np.arange(0,len(amap),1); y = np.arange(0,len(amap),1)
xv, yv = np.meshgrid(x, y)
delX = xv - center[0]
delY = yv - center[1]


from scipy import ndimage

# # test
# component_test = 10**v[0] * np.exp((delX**2 + delY**2 / 0.3**2) /-2 / sigma_pixel**2)
# component_test = ndimage.rotate(component_test, -(pa-90), mode = 'constant', reshape=False)
# plt.figure(); plt.imshow(component_test); plt.colorbar(); plt.show()

# # eq.(33) of Cappellari (2020)
starmass_map = np.zeros(Sigma_map.shape)
luminosity_map = np.zeros(Sigma_map.shape)
for v in ellipse:
    sigma_pixel = v[1] / arcsec_per_pixel
    component = v[0] * np.exp((delX**2 + delY**2 / v[2]**2) / (-2) / sigma_pixel**2)
    if v[2] != 1.0:
        component = ndimage.rotate(component, -(pa - 90), mode = 'constant', reshape=False)
    starmass_map += ml_par * component
    luminosity_map += component
    # plt.figure(); plt.imshow(component); plt.colorbar(); plt.show()
    # plt.figure(); plt.imshow(starmass_map); plt.colorbar(); plt.show()

# starmass_map = starmass_map * np.cos(np.deg2rad(inclination))
pc2_per_pixel = ( pc_per_arcsec * arcsec_per_pixel )**2
print(np.log10(np.sum(starmass_map * pc2_per_pixel))) # total stellar mass

# test
# np.sum(10**ellipse[:,0]) * ml_par
# np.max(starmass_map)

# for v in ellipse:
#     sigma_pixel = v[1] / arcsec_per_pixel
#     map_pixel_sigma = np.sqrt((delX**2 + delY**2 / v[2]**2) / sigma_pixel**2)
#     plt.imshow(map_pixel_sigma); plt.colorbar(); plt.show()


starmass_radial, starmass_radial_error, starmass_radial_std = radial(r_flat_pc_noi, boundary_radial, starmass_map) # mass surface density # , map_weight= 1./Sigma_error**2
luminosity_radial, luminosity_radial_error, luminosity_radial_std = radial(r_flat_pc_noi, boundary_radial, luminosity_map) 

plt.figure()
plt.scatter(plot_radial, luminosity_radial, )
plt.xlabel(r'R$_{\rm gal}$ (pc)')
plt.ylabel(r'log(L / L$_\odot$ pc$^{-2}$)')
plt.yscale('log')
plt.ylim(5,4e5)
# plt.savefig('/Users/ericliang/n1387/work_pub/plot/stellar_luminosity.pdf')
plt.show()

plt.figure()
plt.scatter(plot_radial, starmass_radial, )
plt.xlabel(r'R$_{\rm gal}$ (pc)')
plt.ylabel(r'$\log(\Sigma$ /M$_\odot$ pc$^{-2}$)')
plt.yscale('log')
# plt.savefig('/Users/ericliang/n1387/work_pub/plot/stellar_mass.pdf')
plt.show()


'''
header = 'r_gal (pc); stellar mass surface density (Msun / pc^2)'
data_this = np.array([plot_radial, starmass_radial]).T
np.savetxt('/Users/liangf/work/measurements_publication/profile_star.txt', data_this, header=header)
'''


''' gas fraction plot '''

Sigma_radial_uc = unumpy.uarray(Sigma_radial, Sigma_radial_error)
Sigma_radial_uc_extra = unumpy.uarray(Sigma_radial, np.sqrt(Sigma_radial_error**2 + (0.3*Sigma_radial)**2))
starmass_radial_uc = unumpy.uarray(starmass_radial, starmass_radial_error)
f_gas = Sigma_radial_uc / (Sigma_radial_uc + starmass_radial_uc)
f_gas_extra = Sigma_radial_uc_extra / (Sigma_radial_uc_extra + starmass_radial_uc)


plt.figure(figsize=(4*1.5,3*1.5))
plt.fill_between(ang2pctrans_inv(plot_radial[:-4]), Q[:-4]+Q_err_new[:-4], Q[:-4]-Q_err_new[:-4], alpha=0.3, color='gray')
plt.errorbar(ang2pctrans_inv(plot_radial[:-4]), Q[:-4], yerr=Q_err[:-4], capsize=3, label=r'Toomre parameter $Q$') # ecolor='black'
plt.axvline(ang2pctrans_inv(300),ls='--',color='k');plt.axvline(ang2pctrans_inv(600),ls='--',color='k') # ,label='Region boundary'
plt.xlabel(r'$R_{\mathrm{gal}} $ (arcsec)');plt.ylabel(r'$Q$')
ax1 = plt.gca()
ax1.tick_params(axis='x', which='both', top=False )
secax2 = ax1.secondary_xaxis('top', functions=(ang2pctrans, ang2pctrans_inv))
secax2.set_xlabel(r'$R_{\mathrm{gal}} $ (pc)')
plt.ylim(0,8); plt.xlim(0,13)
plt.text(ang2pctrans_inv(1050),0.65,'Unstable',color='gray'); plt.text(ang2pctrans_inv(1050),1.15,'Stable',color='gray')
plt.axhline(1,ls='-',color='gray') # ,label=r'$Q$ instability threshold'
plt.errorbar([2000], [5], yerr=[1], capsize=3, color='purple',label=r'Molecular gas fraction $f_{\rm mol}$')
plt.legend()

ax5 = ax1.twinx()
plt.fill_between(ang2pctrans_inv(plot_radial[:-4]), unumpy.nominal_values(f_gas_extra[:-4])+unumpy.std_devs(f_gas_extra[:-4]), unumpy.nominal_values(f_gas[:-4])-unumpy.std_devs(f_gas_extra[:-4]), color='pink', alpha=0.6)
plt.errorbar(ang2pctrans_inv(plot_radial[:-4]), unumpy.nominal_values(f_gas[:-4]), yerr=unumpy.std_devs(f_gas[:-4]),capsize=3,color='purple')
plt.ylabel(r'$f_{\rm mol} \equiv \Sigma_{\rm mol} / (\Sigma_{\rm mol} + \Sigma_\star)$',rotation=270,labelpad=20)
plt.ylim(0,0.05)
# plt.gca().set_yticklabels(['{:.0f}%'.format(x*100) for x in plt.gca().get_yticks()])
# plt.savefig('/Users/ericliang/n1387/work_pub/plot/toomre-fgas.pdf')
plt.show()



# dO_dR = np.diff(Omega) / np.diff(curve_rad)
# A =  - curve_rad_clip / 2. * dO_dR # 
# A_radial = np.interp(plot_radial, curve_rad_clip, A)

# # equivalent
# # A=curve_rad * 0.
# # for i in np.arange(0,curve_rad.shape[0]-1):
# #     A[i] = 0.5 * (curve_v[i]/curve_rad[i] - (curve_v[i+1]-curve_v[i])/(curve_rad[i+1]-curve_rad[i]))

# # A_radial = np.interp(plot_radial, curve_rad, A)

# plt.figure()
# plt.plot(plot_radial,A_radial)
# plt.ylabel(r'$A$ [km/s/pc]')
# plt.xlabel(r'$R_{\mathrm{gal}}$')
# plt.show()



# Compare observation with models

# f=open('/Users/liangf/work/output_publication/NGC1387_gmc_table.csv','r')
# reader = csv.reader(f)
# labels = next(reader, None)
# units = next(reader, None)
# cloud_distance = []; cloud_radius = []; cloud_m = []; cloud_sigma = []; cloud_radius_err = []
# for row in reader:
#     if float(row[4]) > 0:
#         cloud_radius.append(float(row[4]))
#         cloud_radius_err.append(float(row[5]))
#         cloud_sigma.append(float(row[6]))
#         cloud_m.append(float(row[12]))
#         cloud_distance.append(float(row[17]))
# cloud_m = np.array(cloud_m) * 1e5
# cloud_distance = np.array(cloud_distance); cloud_radius =np.array(cloud_radius); cloud_sigma = np.array(cloud_sigma)
# cloud_Sigma = cloud_m / (np.pi * cloud_radius**2)
# cloud_radius_err = np.array(cloud_radius_err)



### new version
hdu = fits.open('/Users/ericliang/n1387/work_pub/measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits')
table = hdu[1].data
resolve = np.array(table['resolve_spatial'] * table['resolve_spectral'], dtype=bool)
cloud_distance = table['r_gal'][resolve]
cloud_radius = table['RADRMS_EXTRAP_DECONV'][resolve]
cloud_radius_err = table['RADRMS_EXTRAP_DECONV_UC'][resolve] * cloud_radius
cloud_m = table['mass_extrap'][resolve]
cloud_sigma = table['vrms_extrap_deconv'][resolve]
density = cloud_m / (np.pi * cloud_radius**2)

import lfh
xpoints, ypoints, xerror, yerror = lfh.bin_average(cloud_distance,np.log10(cloud_m),step=10,mode='median')

### needs update, above


theory_radius = 0.19 * Sigma_radial * g_const2 / A_radial**2 # eq 38 of N404 paper
theory_m = Sigma_radial ** 3 * g_const2**2 / 4. / A_radial**4 # eq 17


plt.figure()
plt.axvline(300,ls='--',color='gray');plt.axvline(600,ls='--',color='gray')
plt.plot(plot_radial, np.log10(theory_m), label=r'$M_{\mathrm{c,coll}}$')
plt.scatter(cloud_distance, np.log10(cloud_m), label=r'$M_{\mathrm{clump}}$',s=1)
plt.errorbar(xpoints, ypoints, yerror, xerror, label=r'$M_{\mathrm{clump}}$')
plt.scatter(xpoints, ypoints, marker='x',color='r',s=60)
plt.xlabel(r'$R_{\mathrm{gal}}$')
plt.ylabel(r'$\log(M_{\mathrm{clump}}/{\rm M}_{\odot})$')
plt.legend()
plt.xlim(0,1000)#;plt.ylim(4.4,6.5)
plt.show()

# M = 4 dex M_sun for most regions with R_gal > 400 pc
# Sigma_c = (10**4 / (np.pi * 3.5**2)) = 260 M_sun / pc^2 # normal Sigma_c in N1387
# I = Sigma_c / alpha = 260 / 5.0 = 52 K km/s
np.median( (g_const2 * Sigma_radial / A_radial)[plot_radial>400] ) # =  (max(sigma))_median = 2.9 km/s # sigma eq 26
np.median( ( (2 * g_const2**2 * Sigma_radial**2 / A_radial * 3.5)**(1/3) )[plot_radial>400] ) #  = 2.1 km/s # sigma eq 26
# Iv = I / (2.355*sigma) = 10.5 K


import lfh
xpoints, ypoints, xerror, yerror = lfh.bin_average(cloud_distance,cloud_radius)

plt.figure()
plt.scatter(cloud_distance, cloud_radius, label=r'$R_{\mathrm{obs}}$',s=3,color='blue')
plt.plot(plot_radial, theory_radius, label=r'$R_{\mathrm{model}}$',color='orange')
plt.errorbar(xpoints, ypoints, yerror, xerror, label=r'$R_{\mathrm{clump}}$')
plt.errorbar(1000,21,yerr=np.median(cloud_radius_err), label='median error',capsize=8,color='red')
plt.axvline(300,ls='--',color='gray');plt.axvline(600,ls='--',color='gray')
plt.axhline(14.083707/2, ls='--',label='Resolution limit',color='gray')
plt.xlabel(r'$R_{\mathrm{gal}}$'); plt.ylabel(r'$R_{\rm c}$ (pc)')
plt.legend()
plt.show()
# R <= 2 pc
print(np.median(theory_radius))



##### cloud model prediction of sigma-R relation

boundary = [0,300,600,1e6] # pc

for i in range(3):
    print(i)
    Sigma_ave = np.nanmean(Sigma_map[  (r_flat_pc > boundary[i] ) & (r_flat_pc < boundary[i+1])])
    A_radial_cut = A_radial[(plot_radial > boundary[i] ) & (plot_radial < boundary[i+1])]
    plot_radial_cut = plot_radial[ (plot_radial > boundary[i] ) & (plot_radial < boundary[i+1]) ]
    A_ave = np.sum(A_radial_cut * plot_radial_cut) / np.sum(plot_radial_cut)

    slope = (2 * Sigma_ave**2 * g_const2**2/A_ave)**(1./3.)
    constant = Sigma_ave * g_const2 / A_ave
    critical = Sigma_ave * g_const2 / (2 *A_ave**2)
    print(r'$(2\Sigma ^2 G^2/A)^{1/3}$ =', slope)
    print(r'$\Sigma G /A$ =', constant )
    print(r'$L_D = \lambda_{coll} = $', critical )

    xdata1 = np.linspace(10**0.5, critical, 1000)
    xdata2 = np.linspace(critical, 100, 1000)
    select = (cloud_distance > boundary[i]) & (cloud_distance < boundary[i+1])
    plt.figure()
    plt.scatter(np.log10(cloud_radius[select]),np.log10(cloud_sigma[select]),  s=1)
    plt.plot(np.log10(xdata1), np.log10(xdata1**(1./3.)*slope),ls='--',color='k')
    plt.plot(np.log10(xdata2), np.log10( np.full( len(xdata2), constant)), ls='--',color='k')
    plt.xlabel('log(R/pc)')
    plt.ylabel(r'log($\sigma/(\mathrm{km~s^{-1})}$)')
    plt.xlim(0.6,2.0)
    plt.ylim(-0.3,1.8)
    plt.show()


# m_l = 
# mge_file = '/Users/liangf/work/data/NGC1387_mge_gauss.csv'
# mge_component = np.loadtxt(mge_file, delimiter=",")

# Sigma_star = np.zeros(shape)
# for v in mge_component:
#     pixel_scale = 0.04 # arcsec per pixel
#     x, y = np.meshgrid(np.linspace(-1 * pixel_scale * shape[0] / 2., pixel_scale * shape[0] / 2.,shape[0]), np.linspace(-1 * pixel_scale * shape[0] / 2., pixel_scale * shape[0] / 2., shape[1]))
#     d = np.sqrt(x*x+y*y)
#     sigma, mu = v[1], 0.0
#     g = v[0] * m_l * solar_mass *  np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )

#     Sigma_star += g


# gal_property = np.genfromtxt('/Users/liangf/work/data/NGC1387_galpar2.dat', comments=';',dtype=str)
# center = [int(gal_property[2]), int(gal_property[1])] # integers, [1st index, 2nd index], [RA, DEC]
# bmaj_value = float(gal_property[3])  # arcsec
# bmin_value = float(gal_property[4]) # arcsec
# cellsize = float(gal_property[6]) # arcsec per spaxel
# rest_freq = float(gal_property[7]) # GHz


# # CO 1-0 data
# bmaj_10 = 2.85029 * u.arcsec; bmin_10 = 1.93448 * u.arcsec
# beam_area_10 = 1.1330900354567983 * (bmaj_10*bmin_10)
# factor_10 = (1.0 * (u.Jy/beam_area_10).to(u.K, equivalencies=u.brightness_temperature(115.2712*u.GHz)) ).value 

# # 2-1 ACA data
# bmaj_aca = 6.41563 * u.arcsec; bmin_aca = 5.25789 * u.arcsec
# beam_area_aca = 1.1330900354567983 * (bmaj_aca*bmin_aca)
# factor_aca = (1.0 * (u.Jy/beam_area_aca).to(u.K, equivalencies=u.brightness_temperature(freq)) ).value 

# co_factor = 4.35 * ( 1.0626*factor_10 / (21.408*factor_aca) )  # M_solar / pc^2 / (K km/s)


# plt.figure()
# plt.hist(Sigma_kelvin[Sigma_kelvin>0],bins='auto',histtype='step')
# plt.xlabel('Integrated brightness temperature [K*km/s]')
# plt.ylabel('Spaxel count')
# plt.show()

# plt.figure();
# plt.hist(Sigma_gas[Sigma_gas>0],bins='auto',histtype='step');
# plt.xlabel(r'Mass density [M$_\odot$ / pc$^2$]')
# plt.ylabel('Spaxel count')
# plt.show()




# def radial(map_quantity, map_distance, array_radial, map_weight=None, method='median'):
#     rad_low = array_radial[0]; rad_high = array_radial[-1]
#     quantity_rad = np.zeros(len(array_radial)-1)
#     error_rad = np.zeros(len(array_radial)-1)
#     quantity_rad_all = [[] for _ in range(len(array_radial)-1)]
#     weight_all = [[] for _ in range(len(array_radial)-1)]
#     for i,v1 in enumerate(map_quantity):
#         for j, v in enumerate(v1):
#             if (map_distance[i,j] >= rad_high) or (map_distance[i,j] <= rad_low):
#                 continue
#             index_this = np.searchsorted(array_radial, map_distance[i,j]) -1
#             quantity_rad_all[index_this].append(v)
#             if map_weight is not None:
#                 weight_all[index_this].append(map_weight[i,j])

#     for i2, v2 in enumerate(quantity_rad_all):
#         v2_array = np.array(v2)
#         weight_this = np.array(weight_all[i2])
#         if map_weight is not None:
#             if method == 'mean':
#                 this = DescrStatsW(v2_array[~np.isnan(v2)], weights=weight_this[~np.isnan(v2)], ddof=0)
#                 quantity_rad[i2] = this.mean
#                 # quantity_rad[i2] = ws.weighted_mean(v2_array[~np.isnan(v2)], weights=weight_this[~np.isnan(v2)]) # tested, equivalent to DescrStatsW
#                 error_rad[i2] = this.std
#                 continue
#             quantity_rad[i2] = ws.weighted_median(v2_array[~np.isnan(v2)], weights=weight_this[~np.isnan(v2)])
#             error_rad[i2] = 1.48260 * ws.weighted_median( np.absolute(v2_array[~np.isnan(v2)] - quantity_rad[i2])   , weights=weight_this[~np.isnan(v2)])
#         else:
#             if method == 'mean':
#                 quantity_rad[i2] = np.nanmean(v2)
#                 error_rad[i2] = np.std(v2_array[~np.isnan(v2)])
#                 continue
#             quantity_rad[i2] = np.nanmedian(v2)
#             error_rad[i2] = median_abs_deviation(v2,nan_policy='omit') * 1.48260

#     return quantity_rad, error_rad


