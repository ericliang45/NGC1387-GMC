import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

from scipy import stats

boundary1 = 300; boundary2 = 600; boundary3 = 1050 # max r_gal of valid clouds 1035 pc
region = [0,boundary1,boundary2,boundary3]

lines = open('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat').readlines()
alpha_co = float(lines[18].strip())

resol_ind = np.loadtxt('/Users/ericliang/n1387/work_pub/measurements/NGC1387-resolved_bool.txt')
resolve_true = resol_ind > 0.5

hdu = fits.open('/Users/ericliang/n1387/work_pub/measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits')
table = hdu[1].data
distance = table['r_gal'][resolve_true]

# radius
radius = table['radrms_extrap_deconv'][resolve_true]

sample1 = radius[distance<boundary1]
sample2 = radius[(distance<boundary2)&(distance>boundary1)]
stats.kstest(sample1, sample2)
# pvalue=0.6843548666220022

sample1 = radius[(distance<boundary2)&(distance>boundary1)]
sample2 = radius[distance>boundary2]
stats.kstest(sample1, sample2)
# pvalue=2.128785107370202e-14

sample1 = radius[distance<boundary1]
sample2 = radius[distance>boundary2]
stats.kstest(sample1, sample2)
# pvalue=6.301132256370416e-11


# mass
mass_raw = table['lum_extrap'] * alpha_co # equivalent to the above, but more natural
mass = mass_raw[resolve_true]

sample1 = mass[distance<boundary1]
sample2 = mass[(distance<boundary2)&(distance>boundary1)]
print(stats.kstest(sample1, sample2))
# pvalue=2.2152796494077106e-07

sample1 = mass[(distance<boundary2)&(distance>boundary1)]
sample2 = mass[distance>boundary2]
print(stats.kstest(sample1, sample2))
# pvalue=3.007241460536983e-22

sample1 = mass[distance<boundary1]
sample2 = mass[distance>boundary2]
print(stats.kstest(sample1, sample2))
# pvalue=1.4298564022597041e-37


# v_sigma
vsigma = table['vrms_extrap_deconv'][resolve_true]

sample1 = vsigma[distance<boundary1]
sample2 = vsigma[(distance<boundary2)&(distance>boundary1)]
print(stats.kstest(sample1, sample2))
# pvalue=1.702548991921994e-20

sample1 = vsigma[(distance<boundary2)&(distance>boundary1)]
sample2 = vsigma[distance>boundary2]
print(stats.kstest(sample1, sample2))
# pvalue=5.148876347572864e-25

sample1 = vsigma[distance<boundary1]
sample2 = vsigma[distance>boundary2]
print(stats.kstest(sample1, sample2))
# pvalue=1.1879492081732439e-51

# mass surface density
density = mass / (np.pi * radius**2)

sample1 = density[distance<boundary1]
sample2 = density[(distance<boundary2)&(distance>boundary1)]
print(stats.kstest(sample1, sample2))
# pvalue=1.3111650681493382e-24

sample1 = density[(distance<boundary2)&(distance>boundary1)]
sample2 = density[distance>boundary2]
print(stats.kstest(sample1, sample2))
# pvalue=4.4502857207815684e-11

sample1 = density[distance<boundary1]
sample2 = density[distance>boundary2]
print(stats.kstest(sample1, sample2))
# pvalue=1.5291138719196553e-32


# alpha_vir

a_vir = table['VIRMASS_EXTRAP_DECONV'][resolve_true] / mass

sample1 = a_vir[distance<boundary1]
sample2 = a_vir[(distance<boundary2)&(distance>boundary1)]
print(stats.kstest(sample1, sample2))
# pvalue=0.030723049652777116

sample1 = a_vir[(distance<boundary2)&(distance>boundary1)]
sample2 = a_vir[distance>boundary2]
print(stats.kstest(sample1, sample2))
# pvalue=3.027349327305773e-07

sample1 = a_vir[distance<boundary1]
sample2 = a_vir[distance>boundary2]
print(stats.kstest(sample1, sample2))
# pvalue=1.2468628207730708e-10


