import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

c_speed = 299792.458 # km/s

''' my cubes '''
# final mask -> Jy km/s 304.15723
# clean mask -> Jy km/s 327.17514

''' cleaning (cubes are all on iMac) '''
# fix: using the clean cube and my own cube
# checked each cleaning actually reached the designed depth

# Benchmark, the cube in use for the paper; cleaned to 1.5 mJy/beam
cube_file = '/Users/liangf/work_pub/data/NGC1387_combine_clean5.image.pbcor.fits'
mask_file = '/Users/liangf/work_pub/data/NGC1387_combine_clean9.mask.fits'

cube, hdr = fits.getdata(cube_file, header=True)
mask, hdr_mask = fits.getdata(mask_file, header=True)
# beam = 1.1331*hdr_mask['bmaj']*hdr_mask['bmin']; 
beam = 1.1331 * 0.167104 * 0.135577
pixel = (hdr['CDELT2'] * 3600)**2 ; dv = -hdr['CDELT3'] / hdr['RESTFRQ'] * c_speed

print( np.nansum(cube * mask) * (pixel / beam) * dv, 'Jy km/s' ) # 327

# Reproducing paper draft value
# cube_file = '/Users/liangf/work_pub/data/cube_flux_depth/reproduce/NGC1387_combine_clean_paper.image.pbcor.fits'
cube_file = '/Users/liangf/work_pub/data/cube_flux_depth-mask5/15/NGC1387_combine_clean_paper.image.pbcor.fits'
mask_file = '/Users/liangf/work_pub/data/NGC1387_combine_clean9.mask.fits'

cube, hdr = fits.getdata(cube_file, header=True)
mask, hdr_mask = fits.getdata(mask_file, header=True)
# beam = 1.1331*hdr_mask['bmaj']*hdr_mask['bmin']; 
beam = 1.1331 * 0.155112 * 0.129264
pixel = (hdr['CDELT2'] * 3600)**2 ; dv = -hdr['CDELT3'] / hdr['RESTFRQ'] * c_speed

print( np.nansum(cube * mask) * (pixel / beam) * dv, 'Jy km/s' ) # 329
# same cleaning mask (No.5), same measurement sets, same depths, same summing mask (No.9)
# NB. have to use the accurate beam size (otherwise there is a difference of 10%)
# Comparing directly the products, eg the (negative) minimum of a cube channel, the maximum of a residual channel (outside clean mask), and the beam, they are very similar, but not exactly the same. Don't know why.
# Other related thoughts:
# Is the stored mask the actual mask used in the last step of cleaning (i.e. the final manual adjustment saved or not)
# Does the interactive determination (i.e. cube in use) make a difference from using the final mask in the non-interactive mode (i.e. reproduced one)

# Those below used mask9 in cleaning
# Cleaned to 1.0 mJy/beam
cube_file = '/Users/liangf/work_pub/data/cube_flux_depth/10/NGC1387_combine_clean10m.image.pbcor.fits'
mask_file = '/Users/liangf/work_pub/data/NGC1387_combine_clean9.mask.fits'

cube, hdr = fits.getdata(cube_file, header=True)
mask, hdr_mask = fits.getdata(mask_file, header=True)
beam = 1.1331*hdr_mask['bmaj']*hdr_mask['bmin']; pixel = (hdr['CDELT2'] * 3600)**2 ; dv = -hdr['CDELT3'] / hdr['RESTFRQ'] * c_speed

print( np.nansum(cube * mask) * (pixel / beam) * dv, 'Jy km/s' ) # 244.3362050538479


# Cleaned to 1.5 mJy/beam
cube_file = '/Users/liangf/work_pub/data/cube_flux_depth/15/NGC1387_combine_clean15m.image.pbcor.fits'
mask_file = '/Users/liangf/work_pub/data/NGC1387_combine_clean9.mask.fits'

cube, hdr = fits.getdata(cube_file, header=True)
mask, hdr_mask = fits.getdata(mask_file, header=True)
beam = 1.1331*hdr_mask['bmaj']*hdr_mask['bmin']; pixel = (hdr['CDELT2'] * 3600)**2 ; dv = -hdr['CDELT3'] / hdr['RESTFRQ'] * c_speed

print( np.nansum(cube * mask) * (pixel / beam) * dv, 'Jy km/s' ) # 281.32789291851526


# Cleaned to 2.0 mJy/beam
cube_file = '/Users/liangf/work_pub/data/cube_flux_depth/20/NGC1387_combine_clean20m.image.pbcor.fits'
mask_file = '/Users/liangf/work_pub/data/NGC1387_combine_clean9.mask.fits'

cube, hdr = fits.getdata(cube_file, header=True)
mask, hdr_mask = fits.getdata(mask_file, header=True)
beam = 1.1331*hdr_mask['bmaj']*hdr_mask['bmin']; pixel = (hdr['CDELT2'] * 3600)**2 ; dv = -hdr['CDELT3'] / hdr['RESTFRQ'] * c_speed

print( np.nansum(cube * mask) * (pixel / beam) * dv, 'Jy km/s' ) # 335.65458419294487


# Cleaned to 3.0 mJy/beam
cube_file = '/Users/liangf/work_pub/data/cube_flux_depth/30/NGC1387_combine_clean30m.image.pbcor.fits'
mask_file = '/Users/liangf/work_pub/data/NGC1387_combine_clean9.mask.fits'

cube, hdr = fits.getdata(cube_file, header=True)
mask, hdr_mask = fits.getdata(mask_file, header=True)
beam = 1.1331*hdr_mask['bmaj']*hdr_mask['bmin']; pixel = (hdr['CDELT2'] * 3600)**2 ; dv = -hdr['CDELT3'] / hdr['RESTFRQ'] * c_speed

print( np.nansum(cube * mask) * (pixel / beam) * dv, 'Jy km/s' ) # 450.65439606372877




''' masking (on Macbook) '''
from scipy.ndimage import binary_erosion as ero
from scipy.ndimage import binary_dilation as dia

cube_file = '/Users/ericliang/n1387/work_pub/data/NGC1387_combine_clean5.image.pbcor.fits'
# cube_file = '/Users/ericliang/n1387/work_pub/data/NGC1387_combine_clean5.image.fits'
beam = 1.1331 * 0.167104 * 0.135577

# cube_file = '/Users/ericliang/n1387/work_pub/data/NGC1387_combine_clean10m.image.pbcor.fits'
# beam = 1.1331 * 0.155112 * 0.129264

# mask_file = '/Users/ericliang/n1387/work_pub/data/NGC1387_combine_clean9.mask.fits'
mask_file = '/Users/ericliang/n1387/work_pub/data/NGC1387_combine_clean5.mask.fits'

cube, hdr = fits.getdata(cube_file, header=True)
mask = fits.getdata(mask_file)
length = cube.shape[0]
pixel = (hdr['CDELT2'] * 3600)**2 ; dv = -hdr['CDELT3'] / hdr['RESTFRQ'] * c_speed
origin = np.sum(mask)
# hdr_mask = fits.getheader('/Users/ericliang/n1387/work_pub/data/NGC1387_combine_clean9.mask.fits')
# beam = 1.1331*hdr_mask['bmaj']*hdr_mask['bmin'];

print( np.nansum(cube * mask) * (pixel / beam) * dv, 'Jy km/s' )
print( np.nansum(cube) * (pixel / beam) * dv, 'Jy km/s' ) # 184.3, significant less than the above. Probably negative bowls?

noise_head = cube[:5,:,:]; noise_tail = cube[-5:,:,:]
noise_both = np.concatenate((noise_head,noise_tail))
noise_cube = np.repeat(noise_both, 10, axis=0)

# test noise statistics
# sigma=np.nanstd(noise_cube) # 0.002127
# np.sum(np.isfinite(noise_cube)) # 55006800, 97.8% of total_pix
# np.sqrt(55006800/16) * sigma * (pixel / beam) * dv # sigma of total flux in Jy km/s = 0.5
# np.sqrt(54986090/16) * sigma * (pixel / beam) * dv # sigma of total flux in Jy km/s = 0.5
# np.nansum(noise_cube[10:90]) # -157
# Looking at per channel:
# np.sum(np.isfinite(noise_cube[0])) # 549805
# np.sqrt(549805/16) * sigma  # 0.39 Jy/beam, sigma of sum(each channel)
# np.nansum(noise_cube[0]) # 4.2, -0.6, -1.2, -10.7, 11.8, -9.3, -9.1, 0.9, 2.6, 8.6; actual sigma is 7.2 Jy/beam, 18 times the ideal prediction

# show dialiation performance
# slice_mask = np.array(mask[50],dtype=bool)
# plt.figure()
# plt.imshow(slice_mask)
# plt.savefig('/Users/ericliang/Desktop/3.pdf')
# plt.show()
# for i in range(40):
#     dia(slice_mask, output=slice_mask)

mask_dia = mask.copy()
mask_ero = mask.copy()
total_pix = np.prod(cube.shape)
big_pix=[]; big_flux=[]
small_pix=[]; small_flux=[]
mask_spec = np.sum(mask,axis=(1,2))
ind=np.where(mask_spec>0)[0]


ind = [9]
total_pix = np.prod(cube.shape) / length
origin = np.sum(mask[ind[0]])


def area2width(area):
    shape_ori = area.shape
    area = area.flatten()
    width = area.copy()
    area = area*total_pix
    for i,v in enumerate(area):
        if v==origin: width[i]=0
        if v>origin: width[i] = np.clip( np.argmin(np.abs(v - big_pix)), a_min=0.1, a_max=len(big_pix) )
        if v<origin: width[i] = np.clip( np.argmin(np.abs(v - small_pix)), a_min=0.1, a_max=len(small_pix)) * (-1)
    return np.reshape(width, shape_ori)

def area2width_inv(width):
    shape_ori = width.shape
    width = width.flatten()
    area = width.copy()
    for i,v in enumerate(width):
        if v>0: area[i] = big_pix[np.clip(int(v),0,len(big_pix)-1)] / total_pix
        if v<0: area[i] = small_pix[np.clip(int(-v),0,len(small_pix)-1)] / total_pix
        if v==0: area[i] = origin / total_pix
    return np.reshape(area, shape_ori)

for i in range(800):
    if i%10==0: print(i)
    for j in ind:
        mask_dia[j] = dia(mask_dia[j])
        if i<5: mask_ero[j] = ero(mask_ero[j])
    # big_pix.append(np.sum(mask_dia))
    # big_flux.append(np.nansum(cube * mask_dia) * (pixel / beam) * dv)
    # big_flux.append(np.nansum(noise_cube * mask_dia) * (pixel / beam) * dv)
    big_pix.append(np.sum(mask_dia[ind[0]]))
    big_flux.append(np.nansum(cube[ind[0]] * mask_dia[ind[0]]) * (pixel / beam) * dv)
    if i<20:
        # small_pix.append(np.sum(mask_ero))
        # small_flux.append(np.nansum(cube * mask_ero) * (pixel / beam) * dv)
        # small_flux.append(np.nansum(noise_cube * mask_ero) * (pixel / beam) * dv)
        small_pix.append(np.sum(mask_ero[ind[0]]))
        small_flux.append(np.nansum(cube[ind[0]] * mask_ero[ind[0]]) * (pixel / beam) * dv)

# plt.figure()
# plt.imshow(mask[50])
# plt.savefig('/Users/ericliang/Desktop/00.pdf')
# plt.show()


plt.figure(figsize=(5,4))
plt.scatter(small_pix/total_pix,small_flux)
plt.scatter(big_pix/total_pix,big_flux)
plt.xlabel('Pixel in mask / total pixel')
plt.ylabel('Integrated flux (Jy km/s)')
plt.text(0.1,0.2,'Line-free Channel #'+str(ind[0]),transform=plt.gca().transAxes, fontsize=15)
ax1 = plt.gca()
ax1.tick_params(axis='x', which='both', top=False )
secax2 = ax1.secondary_xaxis('top', functions=(area2width, area2width_inv))
secax2.set_xlabel('Boundary change (pixel)')
plt.show()

# mask-1.png, dia & ero in 3D. 100 loops (only 5 for erosion). Computation rather slow (~10 min). 
# mask-2.png, dia & ero in 2D. Computation ~2 min. Otherwise same as 1.
# mask-2-2.png, 600 loops, approx loop number shown on top. Otherwise same as 2.
# mask-3.png, replace cube with noise_cube. 300 loops. Otherwise same as 2.
# mask-4.png, applied on deeper-cleaned cube (1 Jy/beam). Otherwise same as 2.
# mask-4-2.png, 400 loops, otherwise same as 4.
# mask-5.png, applied on non-PB-corrected cube. Otherwise same as 2.
# mask-6.png, applied on non-PB-corrected cube. Otherwise same as 3.
# mask-6-2.png, 600 loops. Otherwise same as 6.
# mask-7.png, 400 loops, applied to clean5.pbcor cube, with mask5. -> iterations=45



''' beam & cleaning, run on iMac '''

from scipy.ndimage import binary_dilation as dia

pb_cube = fits.getdata('/Users/liangf/work_pub/data/NGC1387_combine_clean5.pb.fits')
pixel = 0.04**2; dv = 2.
high_beam = 1.1331 * 0.155112 * 0.129264 # from CARTA
low_beam = 1.1331 * 0.411057 * 0.34422

mask_clean, header = fits.getdata('/Users/liangf/work_pub/data/NGC1387_mask_full-wrong.fits',header=True)
# /Users/liangf/work_pub/data/NGC1387_combine_clean5.mask.fits
mask_spec = np.sum(mask_clean,axis=(1,2))
ind=np.where(mask_spec>0)[0]
mask_dia = mask_clean.copy()
for j in ind:
    mask_dia[j] = dia(mask_dia[j], iterations=40)

# # temporary, test new masks
# mask_dia = fits.getdata('/Users/liangf/work_pub/NGC1387_mask_full342.fits')
# # temporary

# mask_dia = np.ones(pb_cube.shape)

# cube_file = '/Users/liangf/work_pub/data/NGC1387_combine_clean5.image.pbcor.fits'
# cube, hdr = fits.getdata(cube_file, header=True)
# print(np.nansum(cube * mask_dia) * (pixel / beam) * dv) # 365.73, consistent

high_res = []; high_cl = []; low_res = []; low_cl = []
for v in ['10','15','20','30']:
    print(v)
    high_res_cube = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5/'+v+'/NGC1387_combine_clean'+v+'m.residual.fits') # Jy/beam
    high_cl_cube = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5/'+v+'/NGC1387_combine_clean'+v+'m.model.fits')  # Jy/pix  
    low_res_cube = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5-low_res/'+v+'/NGC1387_combine_clean'+v+'m.residual.fits')
    low_cl_cube = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5-low_res/'+v+'/NGC1387_combine_clean'+v+'m.model.fits')

    high_res.append( np.nansum(high_res_cube * mask_dia / pb_cube) * (pixel / high_beam) * dv)
    high_cl.append( np.nansum(high_cl_cube * mask_dia / pb_cube) * dv)
    low_res.append(np.nansum(low_res_cube * mask_dia/ pb_cube) * (pixel / low_beam) * dv)
    low_cl.append(np.nansum(low_cl_cube * mask_dia/ pb_cube)  * dv)
high_res = np.array(high_res); high_cl = np.array(high_cl); low_res = np.array(low_res); low_cl = np.array(low_cl); 

# plt.figure()
# plt.scatter(high_res,high_cl,label='High-resolution cubes')
# plt.scatter(low_res,low_cl,label='Low-resolution cubes')
# plt.legend()
# plt.xlim(left=0)
# plt.xlabel('Residual flux (Jy km/s)')
# plt.ylabel('Cleaned flux (Jy km/s)')
# plt.show()

plt.figure()
plt.scatter(high_res,high_cl+high_res,label='High-resolution cubes')
plt.scatter(low_res,low_cl+low_res,label='Low-resolution cubes')
plt.legend()
# plt.xlim(left=0); plt.ylim(bottom=100)
# plt.xlim(right=0); plt.ylim(top=280)
plt.xlabel('Residual flux (Jy km/s)')
plt.ylabel('Cleaned flux + Residual flux (Jy km/s)')
plt.axvline(0,ls='--',color='gray')
# plt.show()
plt.savefig('/Users/liangf/work_pub/temp-total_flux.pdf')

hdu = fits.PrimaryHDU(); hdu.data = mask_dia; hdu.header = header
hdu.writeto('mask-.fits', overwrite=True)


''' mask & cleaning & beam (on iMac), for paper '''

from scipy.ndimage import binary_erosion as ero
from scipy.ndimage import binary_dilation as dia

c_speed = 299792.458 # km/s

pb_cube, hdr = fits.getdata('/Users/liangf/work_pub/data/NGC1387_combine_clean5.pb.fits', header=True) # PB cubes are the same within each batch of the same beams, but have tiny differences across the three batches (current high-res, current low-res, actual cube used)
mask = fits.getdata('/Users/liangf/work_pub/data/NGC1387_mask_full.fits') # copied from Macbook, rms_only-202406 version

model_high_15 = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5/15/NGC1387_combine_clean15m.model.fits')
residual_high_15 = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5/15/NGC1387_combine_clean15m.residual.fits')
model_low_15 = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5-low_res/15/NGC1387_combine_clean15m.model.fits')
residual_low_15 = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5-low_res/15/NGC1387_combine_clean15m.residual.fits')

model_high_30 = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5/30/NGC1387_combine_clean30m.model.fits')
residual_high_30 = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5/30/NGC1387_combine_clean30m.residual.fits')
model_low_30 = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5-low_res/30/NGC1387_combine_clean30m.model.fits')
residual_low_30 = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5-low_res/30/NGC1387_combine_clean30m.residual.fits')


# copied from CARTA
beam_low = 1.1331 * 0.411057 * 0.34422
beam_high = 1.1331 * 0.155112 * 0.129264 # This is slightly different from the actual cube's beam of 0.167105" X 0.135577". All imaging parameters are the same though. It could be due to different versions of casa or carta (where the average beam is extracted). The latest test used casa 6.4.4 and carta 3.0.0 while the actual cube's versions used in late 2020 are unknown.

# pixel = (hdr['CDELT2'] * 3600)**2 ; dv = -hdr['CDELT3'] / hdr['RESTFRQ'] * c_speed
pixel = 0.04**2; dv = 2.0 # better numerical precision


mask_spec = np.sum(mask,axis=(1,2))
ind=np.where(mask_spec>0)[0] # 84 channels
# total_pix = np.sum(np.isfinite(residual_high_15)) * len(ind) / len(mask_spec) # 98% of FoV is finite (PB>0.2) in each channel
total_pix = np.prod(residual_high_15.shape) * len(ind) / len(mask_spec)


mask_dia = mask.copy()
mask_ero = mask.copy()
big_pix=[]; big_flux_high=[]; big_flux_low=[]
small_pix=[]; small_flux_high=[]; small_flux_low=[]
for i in range(200):
    if i%20==0: print(i)

    big_pix.append(np.sum(mask_dia))

    # define 1 as shallower (3.0 mJy/beam), 2 as deeper (1.5 mJy/beam)

    c1 = np.nansum(model_high_30 * mask_dia / pb_cube) * dv
    c2 = np.nansum(model_high_15 * mask_dia / pb_cube) * dv
    r1 = np.nansum(residual_high_30 * mask_dia / pb_cube) * (pixel / beam_high) * dv
    r2 = np.nansum(residual_high_15 * mask_dia / pb_cube) * (pixel / beam_high) * dv
    true_flux = (c1 * r2 - r1 * c2) / (r2 - r1)
    big_flux_high.append(true_flux)

    c1 = np.nansum(model_low_30 * mask_dia / pb_cube) * dv
    c2 = np.nansum(model_low_15 * mask_dia / pb_cube) * dv
    r1 = np.nansum(residual_low_30 * mask_dia / pb_cube) * (pixel / beam_low) * dv
    r2 = np.nansum(residual_low_15 * mask_dia / pb_cube) * (pixel / beam_low) * dv
    true_flux = (c1 * r2 - r1 * c2) / (r2 - r1)
    big_flux_low.append(true_flux)
    
    # if i<10:
    #     small_pix.append(np.sum(mask_ero))

    #     c1 = np.nansum(model_high_30 * mask_ero / pb_cube) * dv
    #     c2 = np.nansum(model_high_15 * mask_ero / pb_cube) * dv
    #     r1 = np.nansum(residual_high_30 * mask_ero / pb_cube) * (pixel / beam_high) * dv
    #     r2 = np.nansum(residual_high_15 * mask_ero / pb_cube) * (pixel / beam_high) * dv
    #     true_flux = (c1 * r2 - r1 * c2) / (r2 - r1)
    #     small_flux_high.append(true_flux)
        
    #     c1 = np.nansum(model_low_30 * mask_ero / pb_cube) * dv
    #     c2 = np.nansum(model_low_15 * mask_ero / pb_cube) * dv
    #     r1 = np.nansum(residual_low_30 * mask_ero / pb_cube) * (pixel / beam_low) * dv
    #     r2 = np.nansum(residual_low_15 * mask_ero / pb_cube) * (pixel / beam_low) * dv
    #     true_flux = (c1 * r2 - r1 * c2) / (r2 - r1)
    #     small_flux_low.append(true_flux)

    if i==1: break # to be used after the optimal mask is determined
    for j in ind: # change the mask after the calculation so that the original mask can yield some values
        mask_dia[j] = dia(mask_dia[j],iterations=40)
        if i<5: mask_ero[j] = ero(mask_ero[j],iterations=3)

mask_final = mask_dia.copy()


plt.figure(figsize=(7,5))

plt.scatter(small_pix/total_pix,small_flux_high,label='High-resolution cube, smaller mask',s=1)
plt.scatter(big_pix/total_pix,big_flux_high,label='High-resolution cube, bigger mask',s=1)

plt.scatter(small_pix/total_pix,small_flux_low,label='Low-resolution cube, smaller mask',s=1)
plt.scatter(big_pix/total_pix,big_flux_low,label='Low-resolution cube, bigger mask',s=1)

plt.xlabel('Pixel in mask / total pixel')
plt.ylabel('Total flux (Jy km/s) after cleaning-depth extrapolation')
plt.ylim(200,300)
plt.legend()
plt.savefig('total_flux-three_factors.pdf')
# plt.show()

# # to save another running, at the interval of 4 pixels; at 7 beam-enlargement, 259.5 Jy km/s; probable range 250-270 Jy km/s
# # sophisticated way: use vertical trend moment to get noise and use horizontal trend moment to get true_flux + noise. Using two different masks, the values are 284-27=257 Jy km/s and 265+4=269 Jy km/s, consistent with the range above.

# big_enlargement = np.arange(0,0+4*len(big_flux_high),4)
# plt.figure(figsize=(7,5))
# plt.scatter(big_enlargement, big_flux_high,label='High-resolution cube',s=1)
# plt.xlabel('Enlargement (pixel) or iteration number')
# plt.ylabel('Total flux (Jy km/s) after cleaning-depth extrapolation')
# plt.legend()
# plt.ylim(200,300)
# plt.savefig('Desktop/total_flux-three_factors2.pdf')


# have chosen iterations=28 for the percentage values in the paper. 
# To recover 260 Jy km/s, the cube in use has a mask with pix=40 enlargement. But 28 v 40 has little difference in the values anyway.
print(big_flux_high)
selected_ind = 1
paper_value = big_flux_high[selected_ind]  # 259.52383045289577 Jy km/s

cube_use = fits.getdata('/Users/liangf/work_pub/data/NGC1387_combine_clean5.image.pbcor.fits')
beam_use = 16.044102258194087 * pixel # 1.1331*0.167*0.136
print( np.nansum(cube_use * mask_final) * (pixel / beam_use) * dv * 0.68, 'Jy km/s') # 254.80; different from 253.5 because here the mask can expand outside 20"X20" box while in pymakeplot it's constrained.
# print(np.sum(mask_final[10]))
# plt.figure();plt.imshow(mask_final[10]);plt.savefig('mask_chann.pdf');plt.close()

print(big_pix/total_pix) # 0.0683, location on the curve
print(28*0.04 / np.sqrt(0.167 * 0.136)) # 7.4 beams; beam size from actual cube
print(28*0.04 * 94) # 105.28 pc
print(40*0.04 / np.sqrt(0.167 * 0.136)) # 10.6 beams
print(40*0.04 * 94) # 105.4 pc

# quantification of the differences

print(big_flux_low[selected_ind]/big_flux_high[selected_ind]) # difference caused by beam difference, 0.17%


cube_high_15 = fits.getdata('/Users/liangf/work_pub/data/cube_flux_depth-mask5/15/NGC1387_combine_clean_paper.image.pbcor.fits')

print( np.nansum(cube_high_15 * mask_final) * (pixel / beam_high) * dv, 'Jy km/s' ) # 380.0574
print( np.nansum(cube_high_15 * mask) * (pixel / beam_high) * dv, 'Jy km/s' ) # 286.13276
print(1 - 286.13276 / 380.0574) # difference of 25% due to mask difference (both without cleaning depth correction)

print(380.0574/paper_value) # difference of 46% due to cleaning depth correction
print(paper_value/380.0574) # factor of 0.683, applied to everything relevant to flux

c1 = np.nansum(model_high_30 * mask / pb_cube) * dv
c2 = np.nansum(model_high_15 * mask / pb_cube) * dv
r1 = np.nansum(residual_high_30 * mask / pb_cube) * (pixel / beam_high) * dv
r2 = np.nansum(residual_high_15 * mask / pb_cube) * (pixel / beam_high) * dv
true_flux = (c1 * r2 - r1 * c2) / (r2 - r1) # 226.726
print(1-true_flux / paper_value) # difference of 12.6% due to mask difference (both with cleaning depth correction)
print(paper_value/true_flux-1) # 14.46% growth percentage on the flux mask-growth curve


# one independent correction per channel
# used for paper

raw_flux = []; extrapo_flux = []; epsilon=[]
for j in ind:
    raw_flux.append(np.nansum(cube_high_15[j] * mask_final[j]) * (pixel / beam_high) * dv)
    c1 = np.nansum(model_high_30[j] * mask_final[j] / pb_cube[j]) * dv
    c2 = np.nansum(model_high_15[j] * mask_final[j] / pb_cube[j]) * dv
    r1 = np.nansum(residual_high_30[j] * mask_final[j] / pb_cube[j]) * (pixel / beam_high) * dv
    r2 = np.nansum(residual_high_15[j] * mask_final[j] / pb_cube[j]) * (pixel / beam_high) * dv
    true_flux = (c1 * r2 - r1 * c2) / (r2 - r1)
    extrapo_flux.append(true_flux)
    epsilon.append((c2-c1)/(r1-r2))
raw_flux = np.array(raw_flux); extrapo_flux = np.array(extrapo_flux); epsilon = np.array(epsilon)

print(np.sum(raw_flux)) # 
print(np.sum(extrapo_flux)) # 
print(np.sum(extrapo_flux) / np.sum(raw_flux)) # 
extrapo_flux/raw_flux

plt.figure()
plt.plot(extrapo_flux/raw_flux);
plt.xlabel('channel number');plt.ylabel('Scaling factor per channel');
plt.ylim(0,0.9)
plt.savefig('factor-channel1.pdf')


# common epsilon applied to each channel individually

c1 = np.nansum(model_high_30 * mask_final / pb_cube) * dv
c2 = np.nansum(model_high_15 * mask_final / pb_cube) * dv
r1 = np.nansum(residual_high_30 * mask_final / pb_cube) * (pixel / beam_high) * dv
r2 = np.nansum(residual_high_15 * mask_final / pb_cube) * (pixel / beam_high) * dv
epsilon = (c2-c1)/(r1-r2)  # 0.17

raw_flux = []; extrapo_flux = []
for j in ind:
    raw_flux.append(np.nansum(cube_high_15[j] * mask_final[j]) * (pixel / beam_high) * dv)
    c2 = np.nansum(model_high_15[j] * mask_final[j] / pb_cube[j]) * dv
    r2 = np.nansum(residual_high_15[j] * mask_final[j] / pb_cube[j]) * (pixel / beam_high) * dv
    true_flux = c2 + epsilon*r2
    extrapo_flux.append(true_flux)
raw_flux = np.array(raw_flux); extrapo_flux = np.array(extrapo_flux)

print(np.sum(raw_flux)) # 366.4
print(np.sum(extrapo_flux)) # 256.6
print(np.sum(extrapo_flux) / np.sum(raw_flux)) # 0.70
print(extrapo_flux/raw_flux) # 

plt.figure()
plt.plot(extrapo_flux/raw_flux);
plt.xlabel('channel number');plt.ylabel('Scaling factor per channel');
plt.ylim(0,0.9)
plt.savefig('factor-channel2.pdf')
