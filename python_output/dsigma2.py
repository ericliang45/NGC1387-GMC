import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr

lines = open('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat').readlines()
alpha_co = float(lines[18].strip())

resol_ind = np.loadtxt('/Users/ericliang/n1387/work_pub/measurements/NGC1387-resolved_bool.txt')
resolve = resol_ind > 0.5

hdu_rotation = fits.open('/Users/ericliang/n1387/work_pub/measurements/NGC1387_angmom_comparison.fits')
table_rotation = hdu_rotation[1].data
omega_m = table_rotation['omega_m'][0][resolve]

hdu = fits.open('/Users/ericliang/n1387/work_pub/measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits')
table = hdu[1].data
mass = table['lum_extrap'][resolve] * alpha_co # equiv. to table['mass_extrap'] / table['ALPHA'][0] * alpha_co
distance = table['r_gal'][resolve]
radius = table['RADRMS_EXTRAP_DECONV'][resolve]
id_main = table['peaknum'][resolve]

sigma = table['VRMS_EXTRAP_DECONV'][resolve]
sigma_err = table['VRMS_EXTRAP_DECONV_UC'][resolve] # fractional MAD error
sigma_gs= table['VRMS_GAUSS_EXTRAP_DECONV'][resolve]
sigma_gs_err = table['VRMS_GAUSS_EXTRAP_DECONV_UC'][resolve] # fractional MAD error

diff_obs = sigma**2 - sigma_gs**2
diff_obs_err = np.sqrt( (2*sigma*sigma_err)**2 + (2*sigma_gs*sigma_gs_err)**2 ) 

log_diff_obs_err = diff_obs_err / np.log(10) / diff_obs


pre_table = pd.read_csv('/Users/ericliang/n1387/work_pub/measurements/NGC1387_sig_diff.csv')

prediction = (pre_table['d(sig2)'][1:].to_numpy(dtype=float))[resolve]
prediction_err = (pre_table['  error'][1:].to_numpy(dtype=float))[resolve]

log_prediction_err = prediction_err / np.log(10) / prediction



''' for density_shear '''
''' subsample in Fig 12 '''

density = mass / (np.pi * radius**2)

g_const2 = 4.302e-3 # (pc / Msun) * (km/s)^2
curve_rad, curve_v = np.genfromtxt('/Users/ericliang/n1387/work_pub/data/NGC1387_velocity_curve.dat',comments=';',delimiter=',',unpack=True)

# curve_rad_clip = (curve_rad[1:]+curve_rad[:-1]) / 2.
# Omega = curve_v / curve_rad # km / s / pc

# dO_dR = np.diff(Omega) / np.diff(curve_rad)
# A =  - curve_rad_clip / 2. * dO_dR # 
# A_cloud = np.interp(distance,curve_rad_clip, A  )

# r_t = (g_const2/(2*A_cloud**2))**(1./3.) * mass**(1./3.)
# r_t_err = (1./3.) * (mass_err/mass) * r_t
# r_ratio = radius/r_t
# r_ratio_err = r_ratio * np.sqrt( (radius_err/radius)**2+(r_t_err/r_t)**2 )

# r_shear = 3*np.pi*g_const2*density/(4*A_cloud**2)

''' alternative numerical routine of A_cloud '''
''' conclusion: two methods are in good agreement '''

gal_grid = np.linspace(10,1000,2000)
curve_v_inter = np.interp(gal_grid, curve_rad, curve_v)
gal_grid_clip = (gal_grid[1:]+gal_grid[:-1]) / 2.
Omega_inter = curve_v_inter / gal_grid # km / s / pc
dO_dR_inter = np.diff(Omega_inter) / np.diff(gal_grid)
A_inter =  - gal_grid_clip / 2. * dO_dR_inter # 

A_cloud_2 = np.interp(distance, gal_grid_clip, A_inter  )

# plt.figure()
# plt.scatter(A_cloud,A_cloud_2,s=3)
# plt.plot([0,5],[0,5],ls='--',color='k')
# plt.xlabel('Raw circular velocity array'); plt.ylabel('From interpolated circular velocity array'); 
# plt.gca().set_aspect('equal')
# plt.show()

density_shear = 4 * A_cloud_2**2 * radius / (3 * np.pi * g_const2)



# plt.figure()
# plt.hist(np.log10(density/density_shear)[distance<300],bins='auto',histtype='step',label='Inner') # , cumulative=True
# plt.hist(np.log10(density/density_shear)[distance>300],bins='auto',histtype='step',label='Other') # , cumulative=True
# plt.xlabel(r'$\log(\Sigma_{\rm gas}/\Sigma_{\rm shear})$')
# plt.legend()
# plt.show()

# threshold = 0.6
# np.sum( ((density/density_shear)[distance<300] <1) & ((density/density_shear)[distance<300]>threshold)) # 32
# np.sum( ((density/density_shear)[distance>300] <1) & ((density/density_shear)[distance>300]>threshold)) # 204



''' whole population '''
''' Three colour-coding: c, label, file name '''

a = np.log10(prediction); b =np.log10(diff_obs)
valid = np.isfinite(a) & np.isfinite(b)
a = a[valid]; b = b[valid]
pearsonr(a,b) # PearsonRResult(statistic=0.31779469447585823, pvalue=2.675326301812812e-25)


plt.figure(figsize=(5,3)) # 
plt.scatter(np.log10(prediction),np.log10(diff_obs),s=2,c=distance) # distance  np.log10(mass)  radius
cb = plt.colorbar()
cb.set_label(r'$R_{\rm gal}$ (pc)',rotation=270,labelpad=15) #  '$R_{\rm gal}$ (pc)'  'log$(M_{\rm c}/$M$_\odot$)' '$R_{\rm c}$ (pc)'
plt.errorbar(1,-1,xerr=np.nanmedian(log_prediction_err),yerr=np.nanmedian(log_diff_obs_err))
plt.plot([-5,5],[-5,5],ls='--',color='k')
plt.xlabel(r'log(($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los})_{\rm mod}$ / km$^2$ s$^{-2}$)')
plt.ylabel(r'log($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los}  $ / km$^2$ s$^{-2}$)')
plt.xlim(-4,3);plt.ylim(-4,3)
plt.gca().set_aspect('equal')
plt.savefig('/Users/ericliang/n1387/work_pub/plot/sigma_compare-distance.pdf') # -mass -radius -distance
# plt.show()

plt.figure(figsize=(5,3)) # 
plt.scatter(np.log10(prediction),np.log10(diff_obs),s=2,c=np.log10(mass),vmax=6.2)
cb = plt.colorbar()
cb.set_label(r'$\log(M_{\rm gas}/{\rm M}_\odot)$',rotation=270,labelpad=15)
plt.errorbar(1,-1,xerr=np.nanmedian(log_prediction_err),yerr=np.nanmedian(log_diff_obs_err))
plt.plot([-5,5],[-5,5],ls='--',color='k')
plt.xlabel(r'log(($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los})_{\rm mod}$ / km$^2$ s$^{-2}$)')
plt.ylabel(r'log($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los}  $ / km$^2$ s$^{-2}$)')
plt.xlim(-4,3);plt.ylim(-4,3)
plt.gca().set_aspect('equal')
plt.savefig('/Users/ericliang/n1387/work_pub/plot/sigma_compare-mass.pdf')
# plt.show()

plt.figure(figsize=(5,3)) # 
plt.scatter(np.log10(prediction),np.log10(diff_obs),s=2,c=radius,vmax=40)
cb = plt.colorbar()
cb.set_label(r'$R_{\rm c}$ (pc)',rotation=270,labelpad=15)
plt.errorbar(1,-1,xerr=np.nanmedian(log_prediction_err),yerr=np.nanmedian(log_diff_obs_err))
plt.plot([-5,5],[-5,5],ls='--',color='k')
plt.xlabel(r'log(($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los})_{\rm mod}$ / km$^2$ s$^{-2}$)')
plt.ylabel(r'log($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los}  $ / km$^2$ s$^{-2}$)')
plt.xlim(-4,3);plt.ylim(-4,3)
plt.gca().set_aspect('equal')
plt.savefig('/Users/ericliang/n1387/work_pub/plot/sigma_compare-radius.pdf')
# plt.show()



''' Subsamples '''

selection = radius>30
plt.figure(figsize=(4,3)) # 
plt.scatter(np.log10(prediction[selection & (distance > 300)]),np.log10(diff_obs[selection & (distance > 300)]),s=8, color='deepskyblue', label=r'$R_{\rm gal} \geq 300$ pc')
plt.scatter(np.log10(prediction[selection & (distance < 300)&(omega_m<0.7)]),np.log10(diff_obs[selection & (distance < 300)&(omega_m<0.7)]),s=8, color='purple', label=r'$R_{\rm gal} < 300$ pc' )
plt.scatter(np.log10(prediction[selection & (distance < 300)&(omega_m>0.7)]),np.log10(diff_obs[selection & (distance < 300)&(omega_m>0.7)]),s=8, edgecolor='purple', facecolor='None' )
# plt.text(1,0,r'$R_{\rm c}$ > 30 pc')
plt.plot([-5,5],[-5,5],ls='--',color='k')
plt.xlabel(r'log(($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los})_{\rm mod}$ / km$^2$ s$^{-2}$)')
plt.ylabel(r'log($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los}  $ / km$^2$ s$^{-2}$)')
plt.xlim(-2.2,2.5);plt.ylim(-2.2,2.5)
# plt.legend()
plt.gca().set_aspect('equal')

valid = np.isfinite(np.log10(diff_obs) * np.log10(prediction))
rp,pvalue = pearsonr(np.log10(prediction[selection&valid&(distance>300)]),np.log10(diff_obs[selection&valid&(distance>300)]))
plt.text(0,-0.5,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='dodgerblue',fontsize=9,fontweight='bold')
rp, pvalue = pearsonr(np.log10(prediction[selection & (distance < 300)&(omega_m<0.7)]),np.log10(diff_obs[selection & (distance < 300)&(omega_m<0.7)]))
plt.text(0,-1,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='purple',fontsize=9)

plt.savefig('/Users/ericliang/n1387/work_pub/plot/sigma_compare-sub-radius.pdf')
# plt.show()



selection = mass > 1e6
plt.figure(figsize=(5,3)) # 
plt.scatter(np.log10(prediction[selection & (distance > 300)]),np.log10(diff_obs[selection & (distance > 300)]),s=8, color='deepskyblue', label=r'$R_{\rm gal} \geq 300$ pc')
plt.scatter(np.log10(prediction[selection & (distance < 300)]),np.log10(diff_obs[selection & (distance < 300)]),s=8, color='purple', label=r'$R_{\rm gal} < 300$ pc' )
# plt.text(1,0, r'log$(M_{\rm gas}/$M$_\odot$) > 6')
plt.plot([-5,5],[-5,5],ls='--',color='k')
plt.xlabel(r'log(($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los})_{\rm mod}$ / km$^2$ s$^{-2}$)')
plt.ylabel(r'log($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los}  $ / km$^2$ s$^{-2}$)')
plt.xlim(-2.2,2.5);plt.ylim(-2.2,2.5)
# plt.legend(loc=4)
plt.gca().set_aspect('equal')

valid = np.isfinite(np.log10(diff_obs) * np.log10(prediction))
rp,pvalue = pearsonr(np.log10(prediction[selection&valid&(distance>300)]),np.log10(diff_obs[selection&valid&(distance>300)]))
plt.text(0,-0.5,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='dodgerblue',fontsize=9,fontweight='bold')
rp,pvalue = pearsonr(np.log10(prediction[selection & (distance < 300)]),np.log10(diff_obs[selection & (distance < 300)]))
plt.text(0,-1,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='purple',fontsize=9)

plt.savefig('/Users/ericliang/n1387/work_pub/plot/sigma_compare-sub-mass.pdf')
# plt.show()


selection = (density > 0.6*density_shear) & (density < density_shear)
plt.figure(figsize=(5,3)) # 
plt.scatter(np.log10(prediction[selection & (distance > 300)]),np.log10(diff_obs[selection & (distance > 300)]),s=8, color='deepskyblue', label=r'$R_{\rm gal} \geq 300$ pc')
plt.scatter(np.log10(prediction[selection & (distance < 300)]),np.log10(diff_obs[selection & (distance < 300)]),s=8, color='purple', label=r'$R_{\rm gal} < 300$ pc' )
# plt.text(1,0,r'$\Sigma_{\rm c} < \Sigma_{\rm shear}$')
# plt.text(0.7,-0.3,r'$ 0.6\,\Sigma_{\rm shear} < \Sigma_{\rm c} < \Sigma_{\rm shear}$')
plt.plot([-5,5],[-5,5],ls='--',color='k')
plt.xlabel(r'log(($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los})_{\rm mod}$ / km$^2$ s$^{-2}$)')
plt.ylabel(r'log($\sigma^2_{\rm obs,los} - \sigma^2_{\rm gs,los}  $ / km$^2$ s$^{-2}$)')
plt.xlim(-2.2,2.5);plt.ylim(-2.2,2.5)
# plt.legend()
plt.gca().set_aspect('equal')

valid = np.isfinite(np.log10(diff_obs) * np.log10(prediction))
rp,pvalue = pearsonr(np.log10(prediction[selection&valid&(distance>300)]),np.log10(diff_obs[selection&valid&(distance>300)]))
plt.text(0,-0.5,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='dodgerblue',fontsize=9,fontweight='bold')
rp, pvalue = pearsonr(np.log10(prediction[selection & (distance < 300)]),np.log10(diff_obs[selection & (distance < 300)]))
plt.text(0,-1,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='purple',fontsize=9)

plt.savefig('/Users/ericliang/n1387/work_pub/plot/sigma_compare-sub-density_shear.pdf')
# plt.show()
