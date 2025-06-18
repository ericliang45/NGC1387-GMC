import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import subprocess
from scipy.io import readsav
from scipy.stats import mode

muse = readsav('FCC184_MUSE_gas_measurements.sav',python_dict=True,verbose=True)
bpt = muse['bpt_class']
gas = np.array(muse['oabundance_d16'])
gas[~np.isfinite(gas)] = 0
gas[np.isnan(gas)] = 0
gas_o3n2 = np.array(muse['oabundance_o3n2'])
gas_o3n2[~np.isfinite(gas_o3n2)] = 0
gas_o3n2[np.isnan(gas_o3n2)] = 0

ha_corrected = muse['f_ha']
yn_balmer = muse['yn_balmer']
yn_ha = muse['yn_ha']

gas_sf = gas[bpt == 80]
gas_sf_mean = np.nanmean(gas_sf[gas_sf>0]) # 9.11, expected 8.83 from 2022arXiv220204128L, exactly the same calibration
gas_o3n2_sf = gas_o3n2[bpt == 80]
gas_o3n2_sf_mean = np.nanmean(gas_o3n2[gas_o3n2>0]) # 8.76, expected 8.76 from Pettini & Pagel 2004

x = muse['x']
y = muse['y']

# mode(np.diff(x))

# 50" each side, 0.02" per pixel
gas_2d = np.zeros([int(100/0.2),int(100/0.2)])
gas_o3n2_2d = np.zeros([int(100/0.2),int(100/0.2)])
bpt_2d = np.zeros([int(100/0.2),int(100/0.2)])
x_2d = np.zeros([int(100/0.2),int(100/0.2)])
y_2d = np.zeros([int(100/0.2),int(100/0.2)])
ha_corrected_2d = np.zeros([int(100/0.2),int(100/0.2)])
length = len(x)

for i in range(length):
    x_this = int( (x[i] + 50) / 0.2 ) 
    y_this = int( (y[i] + 50) / 0.2 )
    bpt_2d[y_this,x_this] = bpt[i]
    if 1:
      if bpt[i] == 80:
        gas_o3n2_2d[y_this,x_this] = gas_o3n2[i]
        gas_2d[y_this,x_this] = gas[i]
        if yn_balmer[i] == 1 and yn_ha[i] == 1:
            ha_corrected_2d[y_this,x_this] = ha_corrected[i]


plt.figure()
plt.imshow(gas_2d, vmin=8.5, vmax=9.2, extent=[-50,50,-50,50])
plt.xlim(-10,10);plt.ylim(-10,10)
cb = plt.colorbar()
cb.set_label('12 + log(O/H) (Dopita+16)')
plt.title('Median metallicity (Dopita+16): %.2f' %np.median(gas_2d[gas_2d>0]))
plt.xlabel('RA offset (arcsec)')
plt.ylabel('Dec offset (arcsec)')
plt.show()

plt.figure()
plt.imshow(gas_o3n2_2d, vmin=8.5, vmax=9.2, extent=[-50,50,-50,50])
plt.xlim(-10,10);plt.ylim(-10,10)
cb = plt.colorbar()
cb.set_label('12 + log(O/H) (O3N2)')
plt.title('Median metallicity (O3N2): %.2f' %np.median(gas_o3n2_2d[gas_o3n2_2d>0]))
plt.xlabel('RA offset (arcsec)')
plt.ylabel('Dec offset (arcsec)')
plt.show()


plt.figure()
plt.imshow(bpt_2d, extent=[-50,50,-50,50])
plt.xlim(-10,10);plt.ylim(-10,10)
cb = plt.colorbar(ticks=[80,125,210]) # SF = 80, TR=125, Seyfert=180, LINER=210
cb.ax.set_yticklabels(['SF', 'Composite', 'LINER' ])
plt.title('BPT diagram')
plt.xlabel('RA offset (arcsec)')
plt.ylabel('Dec offset (arcsec)')
plt.show()



plt.figure()
plt.imshow(np.log10(ha_corrected_2d), vmin=3.3,vmax=4.3, extent=[-50,50,-50,50])
plt.xlim(-10,10);plt.ylim(-10,10)
cb = plt.colorbar()
cb.set_label(r'log(Corrected Ha flux / 10$^{-20}$ erg s$^{-1}$ cm$^{-2}$)',rotation=270,labelpad=20)
# plt.title('Median metallicity (Dopita+16): %.2f' %np.median(gas_2d[gas_2d>0]))
plt.xlabel('RA offset (arcsec)')
plt.ylabel('Dec offset (arcsec)')
plt.show()