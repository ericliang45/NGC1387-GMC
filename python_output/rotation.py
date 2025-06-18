import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition # inset_axes , mark_inset
from scipy.stats import pearsonr

boundary1 = 300; boundary2 = 600
# dist = 19.3e6

lines = open('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat').readlines()
dist = float(lines[14].strip())
alpha_co = float(lines[18].strip())

hdu = fits.open('/Users/ericliang/n1387/work_pub/measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits')
table = hdu[1].data
# resolve = np.array(table['resolve_spatial'] * table['resolve_spectral'], dtype=bool)
resol_file = '/Users/ericliang/n1387/work_pub/measurements/NGC1387-resolved_bool.txt'
resol_ind = np.loadtxt(resol_file)
resolve = resol_ind > 0.5

mad_to_sig = 1.48
distance = table['r_gal'][resolve]
radius = table['RADRMS_EXTRAP_DECONV'][resolve]
radius_err = table['RADRMS_EXTRAP_DECONV_UC'][resolve] * radius #  * mad_to_sig
mass = table['lum_extrap'][resolve] * alpha_co
mass_err = table['mass_extrap_uc'][resolve] * mass #  * mad_to_sig
linewidth = table['vrms_extrap_deconv'][resolve]
density = mass / (np.pi * radius**2)

Mvir = (table['VIRMASS_EXTRAP_DECONV']) # tested, the same as vrms_extrap_deconv**2 * radrms_extrap_deconv / (b*G); b=0.2, G=4.302e-3
Mlum = (table['mom0_extrap']*table['chanwidth_kms']*(table['pcperpix'])**2*alpha_co) # tested, the same as lum_extrap * alpha_co
alpha_vir = Mvir/Mlum
alpha_vir_err = np.sqrt( table['VIRMASS_EXTRAP_DECONV_UC']**2 + table['lum_extrap_uc']**2 ) * alpha_vir #  * mad_to_sig

log_alpha = np.log10(alpha_vir[resolve])
log_alpha_err = (1./np.log(10.)*alpha_vir_err/alpha_vir)[resolve]


hdu_rotation = fits.open('/Users/ericliang/n1387/work_pub/measurements/NGC1387_angmom_comparison.fits')
table_rotation = hdu_rotation[1].data
omega_o = table_rotation['omega_o'][0][resolve]
theta_o = table_rotation['theta_o'][0][resolve]
omega_m = table_rotation['omega_m'][0][resolve]
theta_m = table_rotation['theta_m'][0][resolve]
err_omega_o = table_rotation['err_omega_o'][0][resolve]
err_theta_o = table_rotation['err_theta_o'][0][resolve]
err_omega_m = table_rotation['err_omega_m'][0][resolve]
err_theta_m = table_rotation['err_theta_m'][0][resolve]


def ang2pctrans(val): # angle in arcsec, dist in pc
    return val*4.84*dist/1e6

def ang2pctrans_inv(val): # angle in arcsec, dist in pc
    return val/(4.84*dist/1e6)



''' r_t, Sigma_shear calculation '''
g_const2 = 4.302e-3 # (pc / Msun) * (km/s)^2
curve_rad, curve_v = np.genfromtxt('/Users/ericliang/n1387/work_pub/data/NGC1387_velocity_curve.dat',comments=';',delimiter=',',unpack=True)

gal_grid = np.linspace(10,1000,2000)
curve_v_inter = np.interp(gal_grid, curve_rad, curve_v)
gal_grid_clip = (gal_grid[1:]+gal_grid[:-1]) / 2.
Omega_inter = curve_v_inter / gal_grid # km / s / pc
dO_dR_inter = np.diff(Omega_inter) / np.diff(gal_grid)
A_inter =  - gal_grid_clip / 2. * dO_dR_inter # 
A_cloud = np.interp(distance, gal_grid_clip, A_inter  )


# curve_rad_clip = (curve_rad[1:]+curve_rad[:-1]) / 2.
# Omega = curve_v / curve_rad # km / s / pc
# dO_dR = np.diff(Omega) / np.diff(curve_rad)
# A =  - curve_rad_clip / 2. * dO_dR # 
# A_cloud = np.interp(distance,curve_rad_clip, A  )

r_t = (g_const2/(2*A_cloud**2))**(1./3.) * mass**(1./3.)
r_t_err = (1./3.) * (mass_err/mass) * r_t
r_ratio = radius/r_t
r_ratio_err = r_ratio * np.sqrt( (radius_err/radius)**2+(r_t_err/r_t)**2 )

r_shear = 3*np.pi*g_const2*density/(4*A_cloud**2)
density_shear = 4 * A_cloud**2 * radius / (3 * np.pi * g_const2)


''' tidal radius plot '''

plt.figure(figsize=(5,4))
plt.scatter(r_t,radius,c=distance,s=2)
cb = plt.colorbar()
cb.set_label(r'$R_{\rm gal}$ (pc)', rotation=270, labelpad=13)
plt.plot([0,50],[0,50],ls='--',color='gray')
plt.xlim(0,50);plt.ylim(0,50)
plt.xlabel(r'$R_{\rm t}$ (pc)')
plt.ylabel(r'$R_{\rm c}$ (pc)')
plt.gca().set_aspect('equal')
# plt.savefig('/Users/ericliang/n1387/work_pub/plot/tidal_radius.pdf')
plt.show()


plt.figure(figsize=(5,4))
plt.scatter(ang2pctrans_inv(distance), r_ratio,s=2)
# plt.errorbar(distance, r_ratio, yerr=r_ratio_err)
plt.errorbar(ang2pctrans_inv(50), 1, yerr=np.nanmedian(r_ratio_err))
plt.xlim(left=0);plt.ylim(bottom=0)
plt.axhline(1.0,ls='--',color='k')
plt.xlabel(r'$R_{\rm gal}$ (arcsec)')
plt.ylabel(r'$R_{\rm c}/R_{\rm t}$')
plt.axvline(ang2pctrans_inv(300),color='gray');plt.axvline(ang2pctrans_inv(600),color='gray')
plt.text(ang2pctrans_inv(150),6,'Inner',horizontalalignment='center')
plt.text(ang2pctrans_inv(450),6,'Intermediate',horizontalalignment='center')
plt.text(ang2pctrans_inv(800),6,'Outer',horizontalalignment='center')
ax1 = plt.gca()
ax1.tick_params(axis='x', which='both', top=False )
secax2 = ax1.secondary_xaxis('top', functions=(ang2pctrans, ang2pctrans_inv))
secax2.set_xlabel(r'$R_{\mathrm{gal}} $ (pc)')
plt.savefig('/Users/ericliang/n1387/work_pub/plot/tidal_radius_ratio.pdf')
plt.show()


''' Rc / Rt vs alpha_vir '''

plt.figure(figsize=(5,4))
plt.scatter(log_alpha, r_ratio, s=2, c=distance)
cb = plt.colorbar()
cb.set_label(r'$R_{\rm gal}$ (pc)', rotation=270, labelpad=13)
plt.xlabel(r'$\log(\alpha_{\rm obs,vir})$')
plt.ylabel(r'$R_{\rm c}/R_{\rm t}$')
plt.savefig('/Users/ericliang/n1387/work_pub/plot/tidal_radius_alpha.pdf')
plt.show()



fig,ax = plt.subplots(1,3, figsize=(12,4))
plt.sca(ax[0])
select = distance < boundary1
plt.scatter(log_alpha[select], r_ratio[select], s=2)
plt.xlabel(r'$\log(\alpha_{\rm obs,vir})$'); plt.ylabel(r'$R_{\rm c}/R_{\rm t}$')
plt.axhline(1.0,ls='--',color='k');plt.axvline(0.0,ls='--',color='k')
plt.xlim(-1.1,0.9); plt.ylim(0,8.2)
plt.text(-0.9,6,'Inner GMCs')
plt.errorbar(-0.6,5,yerr=np.nanmedian(r_ratio_err[select]),xerr=np.nanmedian(log_alpha_err[select]))
plt.sca(ax[1])
select = (distance > boundary1) & (distance < boundary2)
plt.scatter(log_alpha[select], r_ratio[select], s=2)
plt.xlabel(r'$\log(\alpha_{\rm obs,vir})$')#; plt.ylabel(r'$R_{\rm c}/R_{\rm t}$')
plt.gca().yaxis.set_ticklabels([]) # only hide tick labels
plt.axhline(1.0,ls='--',color='k');plt.axvline(0.0,ls='--',color='k')
plt.xlim(-1.1,0.9); plt.ylim(0,8.2)
plt.text(-0.9,6,'Intermediate GMCs')
plt.errorbar(-0.6,5,yerr=np.nanmedian(r_ratio_err[select]),xerr=np.nanmedian(log_alpha_err[select]))
plt.sca(ax[2])
select = distance > boundary2
plt.scatter(log_alpha[select], r_ratio[select], s=2)
plt.xlabel(r'$\log(\alpha_{\rm obs,vir})$')#; plt.ylabel(r'$R_{\rm c}/R_{\rm t}$')
plt.gca().yaxis.set_ticklabels([]) # only hide tick labels
plt.xlim(-1.1,0.9); plt.ylim(0,8.2)
plt.text(-0.9,6,'Outer GMCs')
plt.axhline(1.0,ls='--',color='k');plt.axvline(0.0,ls='--',color='k')
plt.errorbar(-0.1,5,yerr=np.nanmedian(r_ratio_err[select]),xerr=np.nanmedian(log_alpha_err[select]))

plt.subplots_adjust(wspace=0)
plt.savefig('/Users/ericliang/n1387/work_pub/plot/tidal_radius_alpha_separate.pdf')
plt.show()


# curve_rad_clip = (curve_rad[1:]+curve_rad[:-1]) / 2.
# Omega = curve_v / curve_rad # km / s / pc
# dO_dR = np.diff(Omega) / np.diff(curve_rad)
# A =  - curve_rad_clip / 2. * dO_dR # 
# A_cloud = np.interp(distance,curve_rad_clip, A  )
# r_t = (g_const2/(2*A_cloud**2))**(1./3.) * mass**(1./3.)
# plt.scatter(r_t_new, r_t,s=2);plt.plot([0,59],[0,59]);plt.show()

''' 
# rotation curve
plt.figure(figsize=(3,3))
ax1 = plt.gca()
plt.plot(curve_rad,Omega, label='$\Omega$')
plt.xlim(0,1000)
plt.ylim(0.1,3)
# plt.yscale('log')
plt.plot(gal_grid_clip,A_inter, label='$A$')
plt.ylabel(r'$\Omega$, $A$ (km s$^{-1}$ pc$^{-1}$)')
plt.xlabel(r'$R_{\rm gal} (pc)$')
plt.legend()
# ax5 = ax1.twinx()
# plt.ylim(0.1,3)
plt.savefig('/Users/ericliang/n1387/work_pub/plot/curve_rotation.pdf')
plt.show()
'''


select = np.ones(len(distance), dtype=bool)
inner = (distance < boundary1) & select
intermediate = (distance > boundary1) & (distance < boundary2)  & select
outer = (distance > boundary2)  & select


''' prograde fractions '''

total = len(theta_m); p_num = np.sum(np.abs(theta_m-theta_o)<90)
print('%d clouds, %d prograde, %.2f' %(total,p_num,p_num/total))

total = len(theta_m[inner]); p_num = np.sum(np.abs(theta_m[inner]-theta_o[inner])<90)
print('%d inner clouds, %d prograde, %.2f' %(total,p_num,p_num/total))

total = len(theta_m[intermediate]); p_num = np.sum(np.abs(theta_m[intermediate]-theta_o[intermediate])<90)
print('%d intermediate clouds, %d prograde, %.2f' %(total,p_num,p_num/total))

total = len(theta_m[outer]); p_num = np.sum(np.abs(theta_m[outer]-theta_o[outer])<90)
print('%d outer clouds, %d prograde, %.2f' %(total,p_num,p_num/total))

# 1079 clouds, 675 prograde, 0.63
# 195 inner clouds, 141 prograde, 0.72
# 448 intermediate clouds, 278 prograde, 0.62
# 436 outer clouds, 256 prograde, 0.59



''' Whole-sample plot '''

pearsonr(omega_m, omega_o) # PearsonRResult(statistic=0.10485184756308286, pvalue=0.0005611926)
pearsonr(omega_m[inner], omega_o[inner]) # PearsonRResult(statistic=0.08484000458884111, pvalue=0.23830396)
pearsonr(omega_m[intermediate], omega_o[intermediate]) # PearsonRResult(statistic=0.24541281146121108, pvalue=1.4361822e-07)
pearsonr(omega_m[outer], omega_o[outer]) # PearsonRResult(statistic=0.01819712252232364, pvalue=0.70475644)


pearsonr(theta_m, theta_o) # PearsonRResult(statistic=0.1317197204959175, pvalue=1.4205705e-05)
pearsonr(theta_m[inner], theta_o[inner]) # PearsonRResult(statistic=0.2991456314980496, pvalue=2.155204e-05)
pearsonr(theta_m[intermediate], theta_o[intermediate]) # PearsonRResult(statistic=0.09412213981426021, pvalue=0.046475783)
pearsonr(theta_m[outer], theta_o[outer]) # PearsonRResult(statistic=0.12884021526886613, pvalue=0.007064776)


select = np.ones(len(distance), dtype=bool)
# select = mass > 1e6 #(radius > 30)  & (distance < 250) # & (omega_m < 0.75)

# select_temp = (radius > 30) & (distance < 250) & (distance > 100)
# print(np.corrcoef(omega_m[select_temp], omega_o[select_temp])[0,1])
# print(np.corrcoef(theta_m[select_temp], theta_o[select_temp])[0,1])
# plt.figure()
# plt.scatter(omega_m[select_temp], omega_o[select_temp], c = distance[select_temp])
# plt.colorbar()
# plt.show()

fig = plt.figure(figsize=(10,4))

ax1 = fig.add_subplot(321)
plt.scatter(omega_m[inner], omega_o[inner], s=1, color='blue')
# plt.errorbar(0.5, 0.12, xerr = np.nanmedian(err_omega_m[inner]), yerr = np.nanmedian(err_omega_o[inner]),  color='k', capsize=0)
plt.plot([0,1.5],[0,1.5],ls='--',color='k')
plt.xlim(0,1.5); plt.ylim(0,0.5)
ax1.set_aspect('equal')

ax2 = fig.add_subplot(323, sharex=ax1)
plt.scatter(omega_m[intermediate], omega_o[intermediate], s=1, color='green')
plt.plot([0,1.5],[0,1.5],ls='--',color='k') # ,label='Line of equality'
# plt.errorbar(0.5, 0.12, xerr = np.nanmedian(err_omega_m[intermediate]), yerr = np.nanmedian(err_omega_o[intermediate]), color='k', capsize=1)
ax2.set_aspect('equal')
plt.xlim(0,1.5); plt.ylim(0,0.5)
plt.scatter([], [],s=3, color='blue',label='Inner GMCs')
plt.scatter([],[],s=3,color='green',label='Intermediate GMCs')
plt.scatter([],[],s=3,color='red',label='Outer GMCs')
plt.legend()

ax3 = fig.add_subplot(325, sharex=ax2)
plt.scatter(omega_m[outer], omega_o[outer], s=1, color='red')
plt.plot([0,1.5],[0,1.5],ls='--',color='k')
# plt.errorbar(0.5, 0.12, xerr = np.nanmedian(err_omega_m[outer]), yerr = np.nanmedian(err_omega_o[outer]), label='Median error', color='k', capsize=2)
plt.xlim(0,1.5); plt.ylim(0,0.5)
# plt.legend()
ax3.set_aspect('equal')
plt.errorbar(1.4,0.1, xerr=np.median(err_omega_m), yerr=np.median(err_omega_o), capsize=0, capthick=0, ecolor='k', elinewidth=0.3)
plt.xlabel(r'$\omega_{\rm mod}$ (km s$^{-1}$ pc$^{-1}$)')

ax5 = fig.add_subplot(121, frameon=False)
plt.ylabel(r'$\omega_{\rm obs}$ (km s$^{-1}$ pc$^{-1}$)', labelpad=-10)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False); plt.gca().minorticks_off()

ax4 = fig.add_subplot(122)
plt.scatter(theta_m[outer], theta_o[outer], s=2, color='red')
plt.scatter(theta_m[intermediate], theta_o[intermediate], s=2, color='green')
plt.scatter(theta_m[inner], theta_o[inner], s=2, color='blue')
plt.xlim(-180,180);plt.ylim(-180,180)
plt.plot([-180,180],[-180,180],ls='--',color='k')
plt.fill_between([-180,180],[-90,270],[-270,90], color='gray',alpha=0.4)
plt.errorbar(140,-150, xerr=np.median(err_theta_m), yerr=np.median(err_theta_o), capsize=0, capthick=0, ecolor='k', elinewidth=0.3)
ax4.set_aspect('equal')
plt.xlabel(r'$\phi_{\rm mod}$ (degree)') # ^{\circ}
plt.ylabel(r'$\phi_{\rm rot}$ (degree)')
# plt.errorbar(120, -60, xerr = np.nanmedian(err_theta_m), yerr = np.nanmedian(err_theta_o), color='k', capsize=2)

plt.subplots_adjust(hspace=0, wspace=0.1)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

# plt.savefig('/Users/ericliang/n1387/work_pub/plot/rotation.pdf')
plt.show()

# print('inner r',np.corrcoef(omega_m[inner], omega_o[inner])[0,1])
# print('intermediate r',np.corrcoef(omega_o[intermediate], omega_o[intermediate])[0,1])
# print('outer r',np.corrcoef(omega_m[outer], omega_o[outer])[0,1])

# print('outer r',np.corrcoef(theta_m[outer], theta_o[outer])[0,1])
# print('intermediate r',np.corrcoef(theta_m[intermediate], theta_o[intermediate])[0,1])
# print('inner r',np.corrcoef(theta_m[inner], theta_o[inner])[0,1])




''' Subsample plot(s) '''

''' radius '''

select = (radius > 30)

### outlier justification
# ind = np.where( select & (distance<300)&(omega_m>0.7))
# print(distance[ind]) # 54, 62 pc
# print(radius[ind]) # 31, 32 pc
# print(np.sort(distance[select & (distance<300)]))
# print(np.sort(radius[select & (distance<300)]))
# histo(distance[select & (distance<300)],range=(0,200));plt.show()
### outlier justification

fig,ax = plt.subplots(1,2,figsize=(8,3))

plt.sca(ax[0])
plt.scatter(omega_m[select & (distance>300)], omega_o[select & (distance>300)], s=5, color='deepskyblue', label=r'$R_{\rm gal} \geq 300$ pc' )
plt.scatter(omega_m[select & (distance<300)], omega_o[select & (distance<300)], s=5, color='purple', label=r'$R_{\rm gal} < 300$ pc' ) # , c=distance[select]
plt.scatter(omega_m[select & (distance<300) & (omega_m>0.7)],omega_o[select & (distance<300) & (omega_m>0.7)],color='white',s=0.5)
plt.plot([0,1.5],[0,1.5],ls='--',color='k')
plt.xlim(0,1.2); plt.ylim(0,1.2)
ax[0].set_aspect('equal')
plt.legend()
plt.xlabel(r'$\omega_{\rm mod}$ (km s$^{-1}$ pc$^{-1}$)')
plt.ylabel(r'$\omega_{\rm obs}$ (km s$^{-1}$ pc$^{-1}$)')
plt.text(0.11,0.8,r'$R_{\rm c}$ > 30 pc')
# plt.text(0.2,0.7,r'$R_{\rm gal}$ > 300 pc')

rp, pvalue = pearsonr(omega_m[select & (distance>300)], omega_o[select & (distance>300)])
plt.text(0.05,0.7,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='dodgerblue',fontsize=8,fontweight='bold')
rp, pvalue = pearsonr(omega_m[select & (distance<300)&(omega_m<0.7)], omega_o[select & (distance<300)&(omega_m<0.7)])
plt.text(0.05,0.6,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='purple',fontsize=8)

print(pearsonr(omega_m[select & (distance<300)&(omega_m<7)], omega_o[select & (distance<300)&(omega_m<7)]))
# PearsonRResult(statistic=0.24723028583684936, pvalue=0.23346858)

plt.sca(ax[1])
plt.scatter(theta_m[select & (distance>300)], theta_o[select & (distance>300)], s=5,  color='deepskyblue')
plt.scatter(theta_m[select & (distance<300)], theta_o[select & (distance<300)], s=5,  color='purple') # , color='blue' c=distance[select]
plt.scatter(theta_m[select & (distance<300) & (omega_m>0.7)],theta_o[select & (distance<300) & (omega_m>0.7)],color='white',s=0.5)
# plt.colorbar()
plt.plot([-180,180],[-180,180],ls='--',color='k')
plt.fill_between([-180,180],[-90,270],[-270,90], color='gray',alpha=0.4)
plt.xlim(-180,180);plt.ylim(-180,180)
ax[1].set_aspect('equal')
plt.xlabel(r'$\phi_{\rm mod}$ (degree)')
plt.ylabel(r'$\phi_{\rm rot}$ (degree)')

rp, pvalue = pearsonr(theta_m[select & (distance>300)], theta_o[select & (distance>300)])
plt.text(30,-35,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='dodgerblue',fontsize=7.5,bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'),fontweight='bold')
rp, pvalue = pearsonr(theta_m[select & (distance<300)&(omega_m<0.7)], theta_o[select & (distance<300)&(omega_m<0.7)])
plt.text(30,-70,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='purple',fontsize=7.5,bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'))

print(pearsonr(theta_m[select & (distance<300)], theta_o[select & (distance<300)]))
# PearsonRResult(statistic=0.3339741311901889, pvalue=0.10275961)

plt.subplots_adjust(wspace=0.0)

plt.savefig('/Users/ericliang/n1387/work_pub/plot/rotation-radius.pdf')
# plt.show()


# Pearson correlation coefficient
pearsonr(omega_m[select & (distance<300)&(omega_m<0.7)], omega_o[select & (distance<300) & (omega_m<0.7)])
pearsonr(theta_m[select & (distance<300)&(omega_m<0.7)], theta_o[select & (distance<300) & (omega_m<0.7)])

pearsonr(omega_m[select & (distance>300)], omega_o[select & (distance>300)])
pearsonr(theta_m[select & (distance>300)], theta_o[select & (distance>300)])

# select = (mass > 1e6) # & (distance > 300) # & (radius < 30)
# print(np.corrcoef(omega_m[select & (distance>300)], omega_o[select & (distance>300)])[0,1])
# select = (radius > 30) # & (distance > 300) 
# print(np.corrcoef(omega_m[select & (distance>300) & np.isfinite(omega_m) & np.isfinite(omega_o)], omega_o[select & (distance>300) & np.isfinite(omega_m) & np.isfinite(omega_o)])[0,1])
# select = (radius > 30) # & (distance > 300) 
# print(np.corrcoef(omega_m[select & (distance<300) & np.isfinite(omega_m) & np.isfinite(omega_o)], omega_o[select & (distance<300) & np.isfinite(omega_m) & np.isfinite(omega_o)])[0,1])

''' mass '''

select = (mass > 1e6)

print(np.sort(distance[select & (distance<300)]))

fig,ax = plt.subplots(1,2,figsize=(8,3))

plt.sca(ax[0])
plt.scatter(omega_m[select & (distance>300)], omega_o[select & (distance>300)], s=5, color='deepskyblue', label=r'$R_{\rm gal} \geq 300$ pc')
plt.scatter(omega_m[select & (distance<300)], omega_o[select & (distance<300)], s=5, color='purple', label=r'$R_{\rm gal} < 300$ pc' ) # , c=distance[select]
plt.plot([0,1.5],[0,1.5],ls='--',color='k')
plt.xlim(0,0.7); plt.ylim(0,0.7)
ax[0].set_aspect('equal')
plt.legend()
plt.xlabel(r'$\omega_{\rm mod}$ (km s$^{-1}$ pc$^{-1}$)')
plt.ylabel(r'$\omega_{\rm obs}$ (km s$^{-1}$ pc$^{-1}$)')
plt.text(0.05,0.475,r'$M_{\rm gas}$ > 10$^{6}$ M$_\odot$')
# plt.text(0.2,0.7,r'$R_{\rm gal}$ > 300 pc')

rp, pvalue = pearsonr(omega_m[select & (distance>300)], omega_o[select & (distance>300)])
plt.text(0.04,0.42,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='dodgerblue',fontsize=8,fontweight='bold')
rp, pvalue = pearsonr(omega_m[select & (distance<300)], omega_o[select & (distance<300)])
plt.text(0.04,0.37,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='purple',fontsize=8)


plt.sca(ax[1])
plt.scatter(theta_m[select & (distance>300)], theta_o[select & (distance>300)], s=5,  color='deepskyblue')
plt.scatter(theta_m[select & (distance<300)], theta_o[select & (distance<300)], s=5,  color='purple')
c=distance[select]
# plt.colorbar()
plt.plot([-180,180],[-180,180],ls='--',color='k')
plt.fill_between([-180,180],[-90,270],[-270,90], color='gray',alpha=0.4)
plt.xlim(-180,180);plt.ylim(-180,180)
ax[1].set_aspect('equal')
plt.xlabel(r'$\phi_{\rm mod}$ (degree)')
plt.ylabel(r'$\phi_{\rm rot}$ (degree)')

rp, pvalue = pearsonr(theta_m[select & (distance>300)], theta_o[select & (distance>300)])
plt.text(20,-50,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='dodgerblue',fontsize=8,bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'),fontweight='bold')
rp, pvalue = pearsonr(theta_m[select & (distance<300)], theta_o[select & (distance<300)])
plt.text(20,-75,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='purple',fontsize=8,bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'))


plt.subplots_adjust(wspace=0.0)

plt.savefig('/Users/ericliang/n1387/work_pub/plot/rotation-mass.pdf')
# plt.show()


# Pearson correlation coefficient
pearsonr(omega_m[select & (distance<300)], omega_o[select & (distance<300) ])
pearsonr(theta_m[select & (distance<300)], theta_o[select & (distance<300) ])

pearsonr(omega_m[select & (distance>300)], omega_o[select & (distance>300)])
pearsonr(theta_m[select & (distance>300)], theta_o[select & (distance>300)])


''' shear density '''

ratio = density / density_shear

log_upper = 0.0
# log_lower = np.log10(0.6);  # current ratio 0.6 ~ 1.0 # previous ratio = 1 ~ 1.58
log_lower = np.log10(np.nanmin(ratio)) # updated to remove the lower boundary
# plt.figure(figsize=(3,3))
# plt.hist(np.log10(ratio),bins='auto', histtype='step')
# # plt.axvline(0,ls='--'); plt.axvline(0.2, ls='--',color='')
# plt.ylabel('Number of GMCs')
# x = np.linspace(log_lower,log_upper,10)
# plt.fill_between(x,np.zeros(len(x)),np.zeros(len(x)) + 1000, alpha=0.2, color='red')
# plt.ylim(0,130)
# plt.xlabel(r'$\log(\Sigma_{\rm c} / \Sigma_{\rm shear})$')
# plt.savefig('/Users/ericliang/n1387/work_pub/plot/density_ratio.pdf')

# select = (density < density_shear) # used density < 100 or density > 500


# select = (ratio > 10**log_lower) & (ratio < 10**log_upper)
select = ratio < 10**log_upper



fig,ax = plt.subplots(1,2,figsize=(8,3)) # , layout='constrained'

ax1 = ax[0]
plt.sca(ax1)
plt.scatter(omega_m[select & (distance>300)], omega_o[select & (distance>300)], s=5, color='deepskyblue', label=r'$R_{\rm gal} \geq 300$ pc')
plt.scatter(omega_m[select & (distance<300)], omega_o[select & (distance<300)], s=5, color='purple', label=r'$R_{\rm gal} < 300$ pc' ) # , c=distance[select]
plt.plot([0,3],[0,3],ls='--',color='k')
plt.xlim(0,1.75); plt.ylim(0,1.75)
plt.legend()
plt.xlabel(r'$\omega_{\rm mod}$ (km s$^{-1}$ pc$^{-1}$)')
plt.ylabel(r'$\omega_{\rm obs}$ (km s$^{-1}$ pc$^{-1}$)')
ax1.set_aspect('equal')
plt.text(0.15,1.2,r'$\Sigma_{\rm gas} < \Sigma_{\rm shear}$') # ,fontsize=9

rp, pvalue = pearsonr(omega_m[select & (distance>300)], omega_o[select & (distance>300)])
plt.text(0.08,1.0,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='dodgerblue',fontsize=7,fontweight='bold')
rp, pvalue = pearsonr(omega_m[select & (distance<300)], omega_o[select & (distance<300)])
plt.text(0.08,0.9,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='purple',fontsize=7)

ax2 = plt.axes([0,0,1,1])
ip = InsetPosition(ax1, [0.58,0.37,0.35,0.35]) # x, y, width, height
ax2.set_axes_locator(ip)
plt.sca(ax2)
plt.hist(np.log10(ratio),bins='auto', histtype='step')
x = np.linspace(log_lower,log_upper,10)
plt.fill_between(x,np.zeros(len(x)) + 1000, alpha=0.2, color='red',edgecolor='None')
plt.ylim(0,130); plt.xlim(left=log_lower)
plt.xlabel(r'$\log(\Sigma_{\rm gas} / \Sigma_{\rm shear})$', fontsize=8,labelpad=0)
plt.ylabel('No. of GMCs', fontsize=7, backgroundcolor='none', labelpad=0)
plt.xticks(fontsize=6); plt.yticks(fontsize=6)
# ax2.get_yticklabels()[1].set_backgroundcolor('red') # The box behind the text too large, interfering with the tick
# plt.setp(ax2.get_yticklabels(), backgroundcolor="red") # Controlling all tick labels; box behind text too large

# tb = fig.get_tightbbox(fig.canvas.get_renderer())
# fig.set_size_inches(tb.width, tb.height)
# This method is interfered by set_aspect('equal'), which changes the plotting area. fig.set_size_inches method is a workaround, but is not compatible with current layout adjustment mechanism.
coords = ax1.transAxes.inverted().transform(ax2.get_tightbbox())
border = 0.02
w, h = coords[1] - coords[0] + 2*border
ax1.add_patch(plt.Rectangle(coords[0]-border, w,h, fc="white", transform=ax1.transAxes, zorder=2)) 

plt.sca(ax[1])
plt.scatter(theta_m[select & (distance>300)], theta_o[select & (distance>300)], s=5,  color='deepskyblue')
plt.scatter(theta_m[select & (distance<300)], theta_o[select & (distance<300)], s=5,  color='purple')
c=distance[select]
# plt.colorbar()
plt.plot([-180,180],[-180,180],ls='--',color='k')
plt.fill_between([-180,180],[-90,270],[-270,90], color='gray',alpha=0.4)
plt.xlim(-180,180);plt.ylim(-180,180)
ax[1].set_aspect('equal')
plt.xlabel(r'$\phi_{\rm mod}$ (degree)')
plt.ylabel(r'$\phi_{\rm rot}$ (degree)')

rp, pvalue = pearsonr(theta_m[select & (distance>300)], theta_o[select & (distance>300)])
plt.text(45,-10,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='dodgerblue',fontsize=6.5,bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'),fontweight='bold')
rp, pvalue = pearsonr(theta_m[select & (distance<300)], theta_o[select & (distance<300)])
plt.text(45,-35,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='purple',fontsize=6.5,bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'))

plt.subplots_adjust(wspace=0.0)
plt.savefig('/Users/ericliang/n1387/work_pub/plot/rotation-density_shear.pdf')
# plt.show()


# Pearson correlation coefficient
pearsonr(theta_m[select & (distance<300)], theta_o[select & (distance<300) ])

pearsonr(omega_m[select & (distance>300)], omega_o[select & (distance>300)])
pearsonr(theta_m[select & (distance>300)], theta_o[select & (distance>300)])




''' tidal radius '''

select = (radius > r_t) # & (distance > 300) # & (radius < 30)

fig,ax = plt.subplots(1,2,figsize=(8,3))

plt.sca(ax[0])
plt.scatter(omega_m[select & (distance>300)], omega_o[select & (distance>300)], s=5, color='deepskyblue', label=r'$R_{\rm gal} \geq 300$ pc')
plt.scatter(omega_m[select & (distance<300)], omega_o[select & (distance<300)], s=5, color='purple', label=r'$R_{\rm gal} < 300$ pc' ) # , c=distance[select]
plt.plot([0,1.5],[0,1.5],ls='--',color='k')
# plt.xlim(0,1.); plt.ylim(0,1.)
ax[0].set_aspect('equal')
plt.legend()
plt.xlabel(r'$\omega_{\rm mod}$ (km s$^{-1}$ pc$^{-1}$)')
plt.ylabel(r'$\omega_{\rm obs}$ (km s$^{-1}$ pc$^{-1}$)')
plt.text(0.05,0.65,r'$R_{\rm c}$ > $R_{\rm t}$')
# plt.text(0.2,0.7,r'$R_{\rm gal}$ > 300 pc')

rp, pvalue = pearsonr(omega_m[select & (distance>300)], omega_o[select & (distance>300)])
plt.text(0.04,0.5,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='dodgerblue',fontsize=8,fontweight='bold')
rp, pvalue = pearsonr(omega_m[select & (distance<300)], omega_o[select & (distance<300)])
plt.text(0.04,0.4,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='purple',fontsize=8)


plt.sca(ax[1])
plt.scatter(theta_m[select & (distance>300)], theta_o[select & (distance>300)], s=5,  color='deepskyblue')
plt.scatter(theta_m[select & (distance<300)], theta_o[select & (distance<300)], s=5,  color='purple')
c=distance[select]
# plt.colorbar()
plt.plot([-180,180],[-180,180],ls='--',color='k')
plt.fill_between([-180,180],[-90,270],[-270,90], color='gray',alpha=0.4)
plt.xlim(-180,180);plt.ylim(-180,180)
ax[1].set_aspect('equal')

rp, pvalue = pearsonr(theta_m[select & (distance>300)], theta_o[select & (distance>300)])
plt.text(5,10,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='dodgerblue',fontsize=8,fontweight='bold')
rp, pvalue = pearsonr(theta_m[select & (distance<300)], theta_o[select & (distance<300)])
plt.text(5,-10,r'$r_{\rm p} = %.2f ~(p=%.3f)$' %(rp,pvalue), color='purple',fontsize=8)

plt.xlabel(r'$\theta_{\rm mod}$ (degree)')
plt.ylabel(r'$\phi_{\rm rot}$ (degree)')

plt.subplots_adjust(wspace=0.1)

plt.savefig('/Users/ericliang/n1387/work_pub/plot/rotation-tidal.pdf')
plt.show(); plt.close()



# ''' r_shear '''

# select = (radius > r_shear)

# fig,ax = plt.subplots(1,2,figsize=(8,3))

# plt.sca(ax[0])
# plt.scatter(omega_m[select & (distance>300)], omega_o[select & (distance>300)], s=5, color='purple', label=r'$R_{\rm gal}$ > 300 pc' )
# plt.scatter(omega_m[select & (distance<300)], omega_o[select & (distance<300)], s=5, color='deepskyblue', label=r'$R_{\rm gal}$ < 300 pc' ) # , c=distance[select]
# plt.plot([0,1.5],[0,1.5],ls='--',color='k')
# plt.xlim(0,1.); plt.ylim(0,1.)
# ax[0].set_aspect('equal')
# plt.legend()
# plt.xlabel(r'$\omega_{\rm mod}$ (km s$^{-1}$ pc$^{-1}$)')
# plt.ylabel(r'$\omega_{\rm obs}$ (km s$^{-1}$ pc$^{-1}$)')
# plt.text(0.05,0.65,r'All with $R_{\rm c}$ > $R_{\rm shear}$')

# plt.sca(ax[1])
# plt.scatter(theta_m[select & (distance>300)], theta_o[select & (distance>300)], s=5,  color='purple')
# plt.scatter(theta_m[select & (distance<300)], theta_o[select & (distance<300)], s=5,  color='deepskyblue') # , color='blue' c=distance[select]
# # plt.colorbar()
# plt.plot([-180,180],[-180,180],ls='--',color='k')
# plt.fill_between([-180,180],[-90,270],[-270,90], color='gray',alpha=0.4)
# plt.xlim(-180,180);plt.ylim(-180,180)
# ax[1].set_aspect('equal')

# plt.xlabel(r'$\theta_{\rm mod}\ (^{\circ})$')
# plt.ylabel(r'$\phi_{\rm rot}\ (^{\circ})$')

# plt.subplots_adjust(wspace=0.1)

# plt.savefig('/Users/ericliang/n1387/work_pub/plot/rotation-shear.pdf')
# plt.show()



# ''' v_sigma '''
# select = (linewidth > 1.5) & (linewidth < 2.5)

# fig,ax = plt.subplots(1,2,figsize=(8,3))

# plt.sca(ax[0])
# plt.scatter(omega_m[select & (distance>300)], omega_o[select & (distance>300)], s=5, color='purple', label=r'$R_{\rm gal}$ > 300 pc' )
# plt.scatter(omega_m[select & (distance<300)], omega_o[select & (distance<300)], s=5, color='deepskyblue', label=r'$R_{\rm gal}$ < 300 pc' ) # , c=distance[select]
# plt.plot([0,1.5],[0,1.5],ls='--',color='k')
# plt.xlim(0,1.); plt.ylim(0,1.)
# ax[0].set_aspect('equal')
# plt.legend()
# plt.xlabel(r'$\omega_{\rm mod}$ (km s$^{-1}$ pc$^{-1}$)')
# plt.ylabel(r'$\omega_{\rm obs}$ (km s$^{-1}$ pc$^{-1}$)')
# plt.text(0.05,0.65,r'2.5 km s$^{-1}$ > $\sigma_{\rm c}$ > 1.5 km s$^{-1}$')

# plt.sca(ax[1])
# plt.scatter(theta_m[select & (distance>300)], theta_o[select & (distance>300)], s=5,  color='purple')
# plt.scatter(theta_m[select & (distance<300)], theta_o[select & (distance<300)], s=5,  color='deepskyblue') # , color='blue' c=distance[select]
# # plt.colorbar()
# plt.plot([-180,180],[-180,180],ls='--',color='k')
# plt.fill_between([-180,180],[-90,270],[-270,90], color='gray',alpha=0.4)
# plt.xlim(-180,180);plt.ylim(-180,180)
# ax[1].set_aspect('equal')
# plt.xlabel(r'$\theta_{\rm mod}\ (^{\circ})$')
# plt.ylabel(r'$\phi_{\rm rot}\ (^{\circ})$')

# plt.subplots_adjust(wspace=0.1)
# plt.savefig('/Users/ericliang/n1387/work_pub/plot/rotation-dispersion_small.pdf')
# plt.show()



'''

# Testing different versions of rotation measurements


hdu_rotation1 = fits.open('/Users/ericliang/n1387/work_pub/plot/NGC1387_angmom_comparison-v1.fits') # original
table_rotation1 = hdu_rotation1[1].data
omega_o1 = table_rotation1['omega_o'][0][resolve]
theta_o1 = table_rotation1['theta_o'][0][resolve]
omega_m1 = table_rotation1['omega_m'][0][resolve]
theta_m1 = table_rotation1['theta_m'][0][resolve]

hdu_rotation2 = fits.open('/Users/ericliang/n1387/work_pub/plot/NGC1387_angmom_comparison-v2.fits') # no smmothing
table_rotation2 = hdu_rotation2[1].data
omega_o2 = table_rotation2['omega_o'][0][resolve]
theta_o2 = table_rotation2['theta_o'][0][resolve]
omega_m2 = table_rotation2['omega_m'][0][resolve]
theta_m2 = table_rotation2['theta_m'][0][resolve]

hdu_rotation3 = fits.open('/Users/ericliang/n1387/work_pub/plot/NGC1387_angmom_comparison-v3.fits') # smoothing 10 times larger PSF
table_rotation3 = hdu_rotation3[1].data
omega_o3 = table_rotation3['omega_o'][0][resolve]
theta_o3 = table_rotation3['theta_o'][0][resolve]
omega_m3 = table_rotation3['omega_m'][0][resolve]
theta_m3 = table_rotation3['theta_m'][0][resolve]

hdu_rotation4 = fits.open('/Users/ericliang/n1387/work_pub/plot/NGC1387_angmom_comparison-v4.fits') # correct smoothing
table_rotation4 = hdu_rotation4[1].data
omega_o4 = table_rotation4['omega_o'][0][resolve]
theta_o4 = table_rotation4['theta_o'][0][resolve]
omega_m4 = table_rotation4['omega_m'][0][resolve]
theta_m4 = table_rotation4['theta_m'][0][resolve]

hdu_rotation5 = fits.open('/Users/ericliang/n1387/work_pub/plot/NGC1387_angmom_comparison-v5.fits') # new fitting
table_rotation5 = hdu_rotation5[1].data
omega_o5 = table_rotation5['omega_o'][0][resolve]
theta_o5 = table_rotation5['theta_o'][0][resolve]
omega_m5 = table_rotation5['omega_m'][0][resolve]
theta_m5 = table_rotation5['theta_m'][0][resolve]


omega_o5_cprops = table['velgrad'][resolve]
theta_o5_cprops = np.rad2deg(table['velposang'][resolve])



plt.figure()
plt.scatter(omega_m4, omega_m5)
plt.plot([0,1.5],[0,1.5])
plt.xlabel('correct smmothing')
plt.ylabel('new fitting')
plt.show()


plt.figure()
plt.scatter(omega_o4, omega_o5)
plt.plot([0,1.5],[0,1.5])
plt.xlabel('correct smmothing')
plt.ylabel('new fitting')
plt.show()


# exactly the same
plt.figure()
plt.scatter(omega_o5_cprops, omega_o5)
plt.plot([0,1.5],[0,1.5])
plt.xlabel('correct smmothing')
plt.ylabel('new fitting')
plt.show()

plt.figure()
plt.scatter(theta_o5_cprops, theta_o5)
plt.plot([0,1.5],[0,1.5])
plt.xlabel('correct smmothing')
plt.ylabel('new fitting')
plt.show()


hdu_rotation6 = fits.open('/Users/ericliang/n1387/work_pub/plot/NGC1387_angmom_comparison-v6.fits') # new fitting
table_rotation6 = hdu_rotation6[1].data
theta_o6_test = table_rotation6['theta_o'][0][resolve]

theta_cprops_test = np.rad2deg(table['velposang'][resolve])

plt.figure()
plt.scatter(theta_o6_test, theta_cprops_test)
# plt.plot([0,1.5],[0,1.5])
plt.xlabel('theta v6')
plt.ylabel('theta cprops')
plt.show()


'''
