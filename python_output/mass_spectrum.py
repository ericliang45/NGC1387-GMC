import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData
from scipy.optimize import curve_fit
from scipy.stats import median_abs_deviation as mad
import warnings

lines = open('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat').readlines()
alpha_co = float(lines[18].strip())

fwhm = 4.0 # km/s
mad_to_sig = 1.4826
boundary1 = 300; boundary2 = 600; boundary3 = 1050 # max r_gal of valid clouds 1035 pc
region = [0,boundary1,boundary2,boundary3]

number_density_flag = 0 # 0 means absolute number of clouds N, 1 means number density in unit of kpc^-2
area_inner = np.pi * (boundary1/1000.)**2
area_middle = np.pi * (boundary2/1000.)**2 - area_inner
area_total = np.pi * (boundary3/1000.)**2
area_outer = area_total - area_middle  # using 1000 pc as the outmost distance, for the number density case

cube = fits.getdata('/Users/ericliang/n1387/work_pub/measurements/NGC1387_CO21_cube_2kms_correct.fits')
mad_corrected_cube = mad(cube, nan_policy='omit', axis=None)
print('MAD (corrected cube)',mad_corrected_cube,'K') # 1.32 K
rms = mad_corrected_cube * mad_to_sig # K, mad of the PB-corrected cube * 1.4826

hdu = fits.open('/Users/ericliang/n1387/work_pub/measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits')
# hdu = fits.open('/Users/liangf/work/measurements_publication/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits')
table = hdu[1].data
beam = table['PC2PERBEAM'][0] # pc^2 per beam
print('Beam size',beam,'pc^2')
print('Spectral width',fwhm,'km/s')

# resolve_true = np.array(table['resolve_spatial'] * table['resolve_spectral'], dtype=bool)
resol_ind = np.loadtxt('/Users/ericliang/n1387/work_pub/measurements/NGC1387-resolved_bool.txt')
unres_ind = np.loadtxt('/Users/ericliang/n1387/work_pub/measurements/NGC1387-unresolved_bool.txt')
select_all = (resol_ind + unres_ind) > 0.5
resolve_true = resol_ind > 0.5

# mass_raw = table['mass_extrap'] / table['ALPHA'][0] * alpha_co
mass_raw = table['lum_extrap'] * alpha_co # equivalent to the above, but more natural
mass_min = np.nanmin(mass_raw[resolve_true])
index = np.argmin(mass_raw[resolve_true])
mass_min_err = table['mass_extrap_uc'][resolve_true][index] * mass_min
print('M(min)=',mass_min,'Msun', 'error:', mass_min_err)

delta = rms*fwhm*alpha_co*beam
mass_limit = delta * 10 + mass_min
print('Delta M',delta,'Msun')
print('Completeness = Mmin + 10*Delta_M = ',mass_limit,'Msun')

mass = mass_raw[select_all]
distance = table['r_gal'][select_all]
mass_error = table['mass_extrap_uc'][select_all] * mass

# np.sum(~np.isfinite(mass)) # checked, 0


# for scipy
def sci_power(x, gamma, m0):
    y = (x / m0)**(gamma+1)
    return y

def sci_power_trunc(x, gamma, m0, n0):
    y = n0*((x/m0)**(gamma+1)-1)
    return y

def sci_log_power_trunc(x, gamma, m0, n0):
    y = np.log10(n0) + np.log10( (10**x / m0)**(gamma+1)  -1)
    return y

def sci_log_power(x, gamma, m0):
    y = (gamma + 1) * (x - np.log10(m0))
    return y


# power, gamma, M0
all_nontrun = [-2.41722342e+00,  2.33184398e+07] # [1.73341294e-02, 1.13184800e+06]
inner_nontrun = [-2.41734381e+00,  1.13888931e+07] # [5.43456528e-02, 1.17241986e+06]
intermediate_nontrun = [-2.50967956e+00,  1.16217160e+07] # [3.28040996e-02, 6.96289619e+05]
outer_nontrun = [-2.82778032e+00,  3.75983978e+06] # [2.67436801e-02, 1.89635636e+05]

# power_trunc, gamma, M0, N0
all_trun = [-1.82019202e+00,  1.53280529e+06,  1.80234594e+02] # [1.27979040e-02, 3.84792765e+03, 5.42054942e+00]
inner_trun =  [-1.66705472e+00,  1.81652607e+06,  6.85702154e+01] # # [5.44173817e-02, 1.60254325e+04, 1.01033925e+01]
intermediate_trun = [-1.81593495e+00,  1.48694507e+06,  9.47550841e+01] # [2.49011848e-02, 5.24423994e+03, 4.90950075e+00]
outer_trun = [-2.12609468e+00,  8.72832116e+05,  5.01471771e+01] # [1.79481102e-02, 3.20412508e+03, 1.84884644e+00]


warnings.filterwarnings("ignore")

plt.close('all')
plt.figure(figsize=(5.0,4.8))

## all
number, bins = plt.hist(mass, bins=np.sort(mass), density=False, cumulative=-1, color='white')[0:2]
number = np.append(number, 1)
if number_density_flag == 1:
    number = number / area_total
# plt.plot(bins, number, label='All clouds',color='k')
caps, bars = plt.errorbar(bins, number, xerr=mass_error, label='All GMCs',color='k',ecolor='gray')[1:3]
[bar.set_alpha(0.2) for bar in bars]; [cap.set_alpha(0.2) for cap in caps]

popt, pcov = curve_fit(sci_log_power, np.log10(bins[bins>mass_limit][:-9]), np.log10(number[bins>mass_limit][:-9]), p0=all_nontrun)
plt.plot(bins, 10**sci_log_power(np.log10(bins), *popt),ls=':',color='k')
print('all single\n value:', repr(popt))
print(' error:', repr(np.sqrt(np.diag(pcov))))

popt, pcov = curve_fit(sci_log_power_trunc, np.log10(bins[bins>mass_limit][:-9]), np.log10(number[bins>mass_limit][:-9]),p0=all_trun) #  
plt.plot(bins, 10**sci_log_power_trunc(np.log10(bins), *popt),ls='--',color='k')
print('all trun\n value:', repr(popt))
print(' error:', repr(np.sqrt(np.diag(pcov))))

### inner
number, bins = plt.hist(mass[distance<boundary1], bins=np.sort(mass[distance<boundary1]), density=False, cumulative=-1, color='white')[0:2]
number = np.append(number, 1)
if number_density_flag == 1:
    number = number / area_inner
plt.plot(bins, number, label='Inner GMCs',color='blue')

popt, pcov = curve_fit(sci_log_power, np.log10(bins[bins>mass_limit][:-3]), np.log10(number[bins>mass_limit][:-3]), p0=inner_nontrun)
plt.plot(bins, 10**sci_log_power(np.log10(bins), *popt),ls=':',color='blue')
print('inner single\n value:', repr(popt))
print(' error:', repr(np.sqrt(np.diag(pcov))))

popt, pcov = curve_fit(sci_log_power_trunc, np.log10(bins[bins>mass_limit][:-3]), np.log10(number[bins>mass_limit][:-3]), p0=inner_trun)
plt.plot(bins, 10**sci_log_power_trunc(np.log10(bins), *popt),ls='--',color='blue')
print('inner trun\n value:', repr(popt))
print(' error:', repr(np.sqrt(np.diag(pcov))))


### intermediate
number, bins = plt.hist(mass[(distance>boundary1)&(distance<boundary2)], bins=np.sort(mass[(distance>boundary1)&(distance<boundary2)]), density=False, cumulative=-1, color='white')[0:2]
number = np.append(number, 1)
if number_density_flag == 1:
    number = number / area_middle
plt.plot(bins, number, label='Intermediate GMCs',color='green')

popt, pcov = curve_fit(sci_log_power, np.log10(bins[bins>mass_limit][:-3]), np.log10(number[bins>mass_limit][:-3]), p0=intermediate_nontrun)
plt.plot(bins, 10**sci_log_power(np.log10(bins), *popt),ls=':',color='green')
print('intermediate single\n value:', repr(popt))
print(' error:', repr(np.sqrt(np.diag(pcov))))

popt, pcov = curve_fit(sci_log_power_trunc, np.log10(bins[bins>mass_limit][:-3]), np.log10(number[bins>mass_limit][:-3]), p0=intermediate_trun)
plt.plot(bins, 10**sci_log_power_trunc(np.log10(bins), *popt),ls='--',color='green')
print('intermediate trun\n value:', repr(popt))
print(' error:', repr(np.sqrt(np.diag(pcov))))


### outer
number, bins = plt.hist(mass[distance>boundary2], bins=np.sort(mass[distance>boundary2]), density=False, cumulative=-1, color='white')[0:2]
number = np.append(number, 1)
if number_density_flag == 1:
    number = number / area_outer
plt.plot(bins, number, label='Outer GMCs',color='red')

popt, pcov = curve_fit(sci_log_power, np.log10(bins[bins>mass_limit][:-3]), np.log10(number[bins>mass_limit][:-3]), p0=outer_nontrun)
plt.plot(bins, 10**sci_log_power(np.log10(bins), *popt),ls=':',color='red')
print('outer single\n value:', repr(popt))
print(' error:', repr(np.sqrt(np.diag(pcov))))

popt, pcov = curve_fit(sci_log_power_trunc, np.log10(bins[bins>mass_limit][:-3]), np.log10(number[bins>mass_limit][:-3]), p0=outer_trun)
plt.plot(bins, 10**sci_log_power_trunc(np.log10(bins), *popt),ls='--',color='red')
print('outer trun\n value:', repr(popt))
print(' error:', repr(np.sqrt(np.diag(pcov))))

plt.plot([],ls=':',color='k',label='Single power law')
plt.plot([],ls='--',color='k',label='Truncated power law')

plt.axvline(mass_limit,color='gray',ls='--')
plt.text(mass_limit*1.05,12,'Completeness\nlimit',rotation=90,fontsize=12)

if number_density_flag == 1: plt.ylabel(r'$n$ (>$M_{\rm gas}$) (kpc$^{-2}$)',fontsize=15)
else: plt.ylabel(r'$N$ (>$M_{\rm gas}$)',fontsize=15)
plt.xlabel(r'$M_{\rm gas}$ (M$_{\odot}$)',fontsize=15)
order = [0,1,2,4,3]
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],  loc=3, frameon=True)
plt.xlim(left=3e4); plt.ylim(bottom=0.9, top=2000) # left=np.nanmin(mass)
plt.xscale('log'); plt.yscale('log')
if number_density_flag == 1: plt.savefig('/Users/ericliang/n1387/work_pub/plot/mass_spectrum-density.pdf')
else: plt.savefig('/Users/ericliang/n1387/work_pub/plot/mass_spectrum.pdf')
# plt.show()




''' from IDL mspecfit '''

#          888         888         888         888
# *********************************************************
# FIT WHOLE GALAXY
# TRUNCATED:
# N0, M0, gamma       324.37677       1271697.7      -1.6661718
# err_N0, err_M0, err_gamma      99.1023      54772.4     0.100554
# POWER_LAW:
# N0, gamma       5748507.2      -3.1890620
# err_N0, err_gamma      389242.    0.0533761
# *********************************************************

#          178         178         178         178
# *********************************************************
# FIT INNER REGION:
# TRUNCATED:
# N0, M0, gamma       479.47411       1432082.1      -1.1775201
# err_N0, err_M0, err_gamma      38.6926      92045.4     0.105960
# POWER_LAW:
# N0, gamma       5398675.9      -2.8752980
# err_N0, err_gamma      686218.    0.0953469
# *********************************************************

#          387         387         387         387
# *********************************************************
# FIT INTERMEDIATE REGION:
# TRUNCATED:
# N0, M0, gamma       259.61479       1225582.8      -1.4943104
# err_N0, err_M0, err_gamma      86.1897      55917.2     0.133272
# POWER_LAW:
# N0, gamma       4751629.8      -3.0807010
# err_N0, err_gamma      371384.    0.0703552
# *********************************************************

         # 323         323         323         323
# *********************************************************
# FIT OUTER REGION:
# TRUNCATED:
# N0, M0, gamma            -NaN            -NaN            -NaN
# err_N0, err_M0, err_gamma     0.378411  1.18858e+06  0.000494611
# POWER_LAW:
# N0, gamma            -NaN            -NaN
# err_N0, err_gamma          NaN          NaN
# *********************************************************



# trunc/non-trunc, log/linear, LSE/ODR

# def power(coeff, x): # coeff = [m0,gamma]
#     y = (x / coeff[0])**(coeff[1]+1)
#     return y

# def power_trunc(coeff, x): # coeff = [n0, m0, gamma]
#     y = coeff[0]*((x/coeff[1])**(coeff[2]+1)-1)
#     return y

# def log_power_trunc(coeff, x): # coeff = [n0, m0, gamma]
#     y = np.log10(coeff[0]) + np.log10( (10**x / coeff[1])**(coeff[2]+1)  -1)
#     return y

# def log_power(coeff, x): # coeff = [m0,gamma]
#     y = (coeff[1] + 1) * (x - np.log10(coeff[0]))
#     return y


# replaced by scipy.optimize.curve_fit

# data = RealData(np.log10(bins[bins>mass_limit][:-1]), np.log10(number[bins>mass_limit][:-1])); model = Model(log_power_trunc)
# odr = ODR(data, model, [ 3.5892425e+01,  2.2269960e+06, -1.8545997e+00])
# odr.set_job(fit_type=2) # 2 is least chi2, 0 is full orthogonal distance regression
# output = odr.run()
# output.pprint() # [ 2.25522052e+01  2.02015603e+06 -2.05058570e+00]
# yn = log_power_trunc(output.beta, np.log10(bins))
# plt.plot(bins,10**yn,label='log truncated (LSE) ')

# data = RealData(np.log10(bins[bins>mass_limit]), np.log10(number[bins>mass_limit])); model = Model(log_power)
# odr = ODR(data, model, [7014156.2,-2.6168948])
# odr.set_job(fit_type=2) # 2 is least chi2, 0 is full orthogonal distance regression
# output = odr.run()
# output.pprint() #  [ 6.57979095e+06 -2.64669662e+00]
# yn = log_power(output.beta, np.log10(bins))
# plt.plot(bins, 10**yn, label='log non-truncated (LSE)')


# not needed for this

# data = RealData(np.log10(bins[bins>mass_limit]), np.log10(number[bins>mass_limit])); model = Model(log_power)
# odr = ODR(data, model, [7014156.2,-2.6168948]) # 1 is free , ifixb=[1,1,1] 
# odr.set_job(fit_type=0) # 2 is least chi2, 0 is full orthogonal distance regression
# output = odr.run()
# output.pprint() #  [ 4.77951137e+06 -2.89074255e+00]
# yn = log_power(output.beta, np.log10(bins))
# plt.plot(bins, 10**yn, label='log non-truncated (ODR)')


# cross-tested with scipy.optimize.curve_fit

# popt, pcov = curve_fit(sci_power, bins[bins>mass_limit], number[bins>mass_limit], p0=[7014156.2, -2.6168948])
# # [ 7.40206603e+07, -1.87564359e+00] # in good agreement

# popt, pcov = curve_fit(sci_power_trunc, bins[bins>mass_limit], number[bins>mass_limit], p0=[ 3.5892425e+01,  2.2269960e+06, -1.8545997e+00])
# # plt.plot(bins, sci_power_trunc(bins, *popt), label='linear truncated (LSE)')
# # [ 1.96299175e+03,  1.42319325e+06, -1.04148062e+00] # ODR package failed

# popt, pcov = curve_fit(sci_log_power_trunc, np.log10(bins[bins>mass_limit][:-1]), np.log10(number[bins>mass_limit][:-1]), p0=[ 3.5892425e+01,  2.2269960e+06, -1.8545997e+00])
# # [ 2.25535310e+01,  2.02013196e+06, -2.05056393e+00] # in good agreement, also must exclude the largest bin

# popt, pcov = curve_fit(sci_log_power, np.log10(bins[bins>mass_limit]), np.log10(number[bins>mass_limit]), p0=[7014156.2, -2.6168948])
# # [ 6.57979130e+06, -2.64669659e+00] # in good agreement



# not preferred

# data = RealData(bins[bins>mass_limit], number[bins>mass_limit]); model = Model(power)
# odr = ODR(data, model, [7014156.2,-2.6168948]) # 1 is free , ifixb=[1,1,1] 
# odr.set_job(fit_type=0) # 2 is least chi2, 0 is full orthogonal distance regression
# output = odr.run()
# output.pprint() #  [ 7.40210931e+07 -1.87564264e+00]
# yn = power(output.beta, bins)
# plt.plot(bins, yn, '+', label='linear non-truncated (ODR)',ls=':')

# data = RealData(bins[bins>mass_limit], number[bins>mass_limit]); model = Model(power)
# odr = ODR(data, model, [7014156.2,-2.6168948]) # 1 is free , ifixb=[1,1,1] 
# odr.set_job(fit_type=2) # 2 is least chi2, 0 is full orthogonal distance regression
# output = odr.run()
# output.pprint() #  [ 7.40211468e+07 -1.87564252e+00]
# yn = power(output.beta, bins)
# plt.plot(bins, yn, label='linear non-truncated (LSE)',ls='--')



# failed 

# linear truncated LSE
# data = RealData(bins[bins>mass_limit], number[bins>mass_limit]); model = Model(power_trunc)
# odr = ODR(data, model, [ 3.5892425e+01,  2.2269960e+06, -1.8545997e+00])
# odr.set_job(fit_type=2) # 2 is least chi2, 0 is full orthogonal distance regression
# output = odr.run()
# output.pprint() #  wrong

# linear truncated ODR
# data = RealData(bins[bins>mass_limit], number[bins>mass_limit]); model = Model(power_trunc)
# odr = ODR(data, model, [ 3.5892425e+01,  2.2269960e+06, -1.8545997e+00])
# odr.set_job(fit_type=0) # 2 is least chi2, 0 is full orthogonal distance regression
# output = odr.run()
# output.pprint() # wrong

# log truncated ODR
# data = RealData(np.log10(bins[bins>mass_limit][:-1]), np.log10(number[bins>mass_limit][:-1])); model = Model(log_power_trunc)
# odr = ODR(data, model, [ 3.5892425e+01,  2.2269960e+06, -1.8545997e+00])
# odr.set_job(fit_type=0) # 2 is least chi2, 0 is full orthogonal distance regression
# output = odr.run()
# output.pprint() # wrong
