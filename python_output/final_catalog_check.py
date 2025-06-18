''' check criteria: contrast (delta), convexity3D, minarea, minvchan '''
''' distinguish resolved and unresolved clouds '''
''' output two lists of boolean values '''

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

minarea = 16 # spaxel
minvchan = 2 # channel
# prop_file = '/Users/ericliang/n1387/work_pub/measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits' # Macbook Air
# assign_file = '/Users/ericliang/n1387/work_pub/measurements/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits'
# data_file = '/Users/ericliang/n1387/work_pub/measurements/NGC1387_CO21_cube_2kms_correct.fits'
# out_reso = '/Users/ericliang/n1387/work_pub/measurements/NGC1387-resolved_bool.txt'
# out_unre = '/Users/ericliang/n1387/work_pub/measurements/NGC1387-unresolved_bool.txt'
# minconvexity = 0.50
# delta = 2.0 # K

minconvexity = 0.45; delta = 3.0; number = '/appendix/1'
# minconvexity = 0.45; delta = 2.0; number = '/appendix/2'
# minconvexity = 0.45; delta = 1.0; number = '/appendix/3'
# minconvexity = 0.55; delta = 3.0; number = '/appendix/4'
# minconvexity = 0.55; delta = 2.0; number = '/appendix/5'
# minconvexity = 0.55; delta = 1.0; number = '/appendix/6'
# minconvexity = 0.50; delta = 3.0; number = '/appendix/7'
# minconvexity = 0.50; delta = 2.0; number = '/appendix/8' # the chosen one
# minconvexity = 0.50; delta = 1.0; number = '/appendix/9'

prop_file = '/Users/liangf/work_pub/measurements'+number+'/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits' # iMac
assign_file = '/Users/liangf/work_pub/measurements'+number+'/NGC1387_CO21_cube_2kms_assign_clfriendtoo.fits' # iMac
data_file = '/Users/liangf/work_pub/measurements/NGC1387_CO21_cube_2kms_correct.fits' # iMac
out_reso = '/Users/liangf/work_pub/measurements'+number+'/NGC1387-resolved_bool.txt'
out_unre = '/Users/liangf/work_pub/measurements'+number+'/NGC1387-unresolved_bool.txt'

catalog = fits.getdata(prop_file)
assign_cube = fits.getdata(assign_file) # assign_cube.shape = (100, 750, 750)
data_cube = fits.getdata(data_file)

convexity3d = catalog['CONVEXITY3D']

resolve_spectral = catalog['RESOLVE_SPECTRAL']
resolve_spatial = catalog['RESOLVE_SPATIAL']
resolve_flag = resolve_spectral * resolve_spatial

cloud_resolved_list = np.zeros(len(catalog),dtype=bool)
cloud_unresol_list = np.zeros(len(catalog),dtype=bool)

contrast_all = []
width_all = []
area_all = []

for i in range(len(catalog)):

    if i%100 ==0: print(i)

    ind_cloud = np.where(assign_cube == i+1)
    spec_ind = ind_cloud[0]
    spec_width = len(np.unique(spec_ind))

    intensity = data_cube[ind_cloud]
    contrast = np.max(intensity) - np.min(intensity)
    if np.isnan(contrast): print('Error error error')
    contrast_all.append(contrast)

    spatial_ind = np.array([ind_cloud[1],ind_cloud[2]]).T
    area = len(np.unique(spatial_ind,axis=0))
    area_all.append(area)
    width_all.append(spec_width)

    convex_this = convexity3d[i]

    if (area >= minarea) and (spec_width>=minvchan) and (convex_this>=minconvexity) and (contrast>=delta):
        if resolve_flag[i]:
            cloud_resolved_list[i] = True
        else:
            cloud_unresol_list[i] = True
    else: # it's not a cloud any more
        # print('Cloud Number,',catalog[i]['peaknum'])
        # print('(area >= minarea) and (spec_width>=minvchan) and (convex_this>=minconvexity) and (contrast>=delta)')
        # print((area >= minarea), (spec_width>=minvchan), (convex_this>=minconvexity), (contrast>=delta))
        # print('\n')
        pass 

np.savetxt(X=cloud_resolved_list, fname=out_reso)
np.savetxt(X=cloud_unresol_list, fname=out_unre)

contrast_all = np.array(contrast_all); area_all = np.array(area_all); width_all = np.array(width_all)

bins = np.linspace(0,15,31)
plt.figure()
plt.hist(contrast_all,bins=bins)
plt.xlabel('Contrast')
plt.title('Min contrast = '+str(delta))
# plt.show()
plt.savefig('measurements'+number+'/contrast.pdf')

plt.figure()
plt.hist(contrast_all[(convexity3d>minconvexity)&(width_all>minvchan)&(area_all>minarea)],bins=bins)
plt.xlabel('Contrast')
plt.title('Min contrast = '+str(delta))
# plt.show()
plt.savefig('measurements'+number+'/contrast_pure_initial_catalog.pdf')


plt.figure()
plt.hist(convexity3d,bins='auto')
plt.xlabel('Convexity')
plt.title('Min convexity = '+str(minconvexity))
# plt.show()
plt.savefig('measurements'+number+'/convexity.pdf')



''' compare Appendix No.8 and the catalogue in paper '''

hdu = fits.open('/Users/ericliang/ox_local/work_pub/measurements/appendix/8/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits')
table = hdu[1].data
radius_err = table['VIRMASS_EXTRAP_DECONV_UC']
cloudid = table['peaknum']
res = np.loadtxt('/Users/ericliang/ox_local/work_pub/measurements/appendix/8/NGC1387-resolved_bool.txt')

ind = np.where(np.isnan(radius_err[res>0.5]))
print(cloudid[res>0.5][ind]) # 522,  600, 1077, 1394


hdu = fits.open('/Users/ericliang/ox_local/work_pub/measurements/NGC1387_CO21_cube_2kms_props_clfriendtoo.fits')
table = hdu[1].data
radius_err = table['VIRMASS_EXTRAP_DECONV_UC']
cloudid = table['peaknum']
res = np.loadtxt('/Users/ericliang/ox_local/work_pub/measurements/NGC1387-resolved_bool.txt')

ind = np.where(np.isnan(radius_err[res>0.5]))
print(cloudid[res>0.5][ind]) # 1394

''' The difference is probably due to the MCMC nature of the error estimate '''

