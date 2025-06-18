import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


''' HST, Figure 1 '''

from astropy import wcs
from astropy.coordinates import Angle
from scipy import ndimage
from matplotlib.colors import SymLogNorm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from matplotlib.patches import Ellipse

# from matplotlib.colors import LogNorm
# from matplotlib.patches import Ellipse
# from astropy import units as u
# from astropy.coordinates import SkyCoord
# from skimage import data
# from skimage.filters import unsharp_mask
# from cmcrameri import cm


def get_header_coord_arrays(hdr):
    wcs_this=wcs.WCS(hdr)
    print('check square pixel, check ra dec order:')
    print(wcs_this)

    xp = np.arange(hdr['NAXIS1'])
    yp = np.ones(len(xp)) * int(hdr['NAXIS2']/2)
    ra = wcs_this.all_pix2world(xp,yp, 0)[0]

    yp = np.arange(hdr['NAXIS2'])
    xp = np.ones(len(yp)) * int(hdr['NAXIS1']/2)
    dec = wcs_this.all_pix2world(xp,yp, 0)[1]

    try:
        delta = hdr['CD2_2'] * 3600
    except:
        delta = hdr['CDELT2'] * 3600
    cd3= np.nan
    v1 = np.nan

    return ra,dec,v1,delta,cd3

def ang2pctrans_inv(val):
    return val/(4.84*gal_distance)

def ang2pctrans(val):
    return val*4.84*gal_distance

def scalebar(ax,loc='lower left'):
    barlength_pc = 500
    barlength_arc=  barlength_pc/(4.84*gal_distance)
    label = '500 pc'
    asb = AnchoredSizeBar(ax.transData,  barlength_arc,   label,  loc=loc,  pad=0.25, borderpad=0.5, sep=5, frameon=False)
    ax.add_artist(asb)


hst_hdu = fits.open('/Users/ericliang/n1387/work_pub/data/hst_10217_09_acs_wfc_f475w_sci.fits')
hst_image = hst_hdu[1].data

xcoord,ycoord,vcoord,cellsize,dv = get_header_coord_arrays(hst_hdu[1].header)

lines = open('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat').readlines()
ra = Angle(lines[10].strip())
dec = Angle(lines[12].strip())
obj_ra = ra.degree
obj_dec = dec.degree
obj_dec_hst = obj_dec - 0.5/3600 # astrometric correction
imagesize = [12.2,12.2] # half side length in arcsec
posang = float(lines[4].strip())
inc = float(lines[6].strip())
gal_distance = float(lines[14].strip()) / 1e6

wx,=np.where(np.abs((xcoord-obj_ra)*3600 * np.cos(np.deg2rad(obj_dec))) <= imagesize[0]) # edited by Eric
wy,=np.where(np.abs((ycoord-obj_dec_hst)*3600) <= imagesize[1])
spatial_trim=[np.min(wx),np.max(wx),np.min(wy),np.max(wy)]        

image_trim=hst_image[spatial_trim[2]:spatial_trim[3],spatial_trim[0]:spatial_trim[1]]
xcoord_trim=xcoord[spatial_trim[0]:spatial_trim[1]]
ycoord_trim=ycoord[spatial_trim[2]:spatial_trim[3]]
ra_relative = (xcoord_trim-obj_ra)*3600 * np.cos(np.deg2rad(obj_dec))
dec_relative = (ycoord_trim-obj_dec_hst)*3600

# HST spatial resolution 0.1", cell_size = 0.05", 20*0.05=1", 10 PSF
hst_median = ndimage.median_filter(image_trim, size=21)
hst_min = ndimage.minimum_filter(image_trim, size=41)
difference_median = image_trim - hst_median
enhanced = image_trim + 4 * difference_median
mini_result = image_trim - hst_min

cloud_hdu = fits.open('/Users/ericliang/n1387/work_pub/plot/moment_maps/rms_only-202406/NGC1387_mom0.fits')
cloud_map = cloud_hdu[0].data
cloud_map[np.isnan(cloud_map)] = 0.
kernel = Gaussian2DKernel(x_stddev=8) # used 4 in 20230101 version
smooth_cloud = convolve(cloud_map, kernel)

xcoord_cloud,ycoord_cloud,vcoord_cloud,cellsize_cloud,dv_cloud = get_header_coord_arrays(cloud_hdu[0].header)
ra_relative_cloud = (xcoord_cloud - obj_ra)*3600 * np.cos(np.deg2rad(obj_dec))
dec_relative_cloud = (ycoord_cloud - obj_dec)*3600 



''' upper panel '''

plt.figure(figsize=(3.5,3.5))
plt.imshow(difference_median,cmap='Greys_r',norm=SymLogNorm(linthresh=0.3, vmin=-50,vmax=50), extent=[ra_relative[0],ra_relative[-1],dec_relative[0],dec_relative[-1]])
ax1 = plt.gca()
ax1.set_aspect('equal')

ax1.set_xticklabels([])
plt.ylabel(r'Dec. - Dec.$\mathrm{_{centre}}$ (arcsec)',fontsize=15)
ax1.tick_params(axis='y', which='both', right=False, ); ax1.tick_params(axis='x', which='both', top=False )
secax2 = ax1.secondary_xaxis('top', functions=(ang2pctrans, ang2pctrans_inv))
secax2.set_xlabel(r'RA - RA$\mathrm{_{centre}}$ (pc)',fontsize=15)
secax = ax1.secondary_yaxis('right', functions=(ang2pctrans, ang2pctrans_inv))
secax.set_ylabel(r'Dec. - Dec.$\mathrm{_{centre}}$ (pc)',fontsize=15,rotation=270,labelpad=20)


circle = [300,600]
for v in circle:
    radius = ang2pctrans_inv(v) #  / self.cellsize
    circle1 = Ellipse(xy=(0,0), width=radius*2, height=2*radius*np.cos(np.deg2rad(inc)), angle=(360-posang)-90, 
        edgecolor='r', fc='None', lw=1,  zorder=10)
    ax1.add_patch(circle1)

# to calibrate orientation
# ax1.contour(ra_relative_cloud, dec_relative_cloud, cloud_mom1, cmap='jet',levels=[1150,1200,1250,1300,1350,1400], zorder=1,alpha=0.4)

scalebar(ax1,)

plt.text(11,10,'$HST$',fontsize=13)

plt.savefig('/Users/ericliang/n1387/work_pub/plot/hst_upper.pdf')
plt.show()



''' lower panel '''
plt.figure(figsize=(3.5,3.5))
plt.imshow(difference_median,cmap='Greys_r',norm=SymLogNorm(linthresh=0.3, vmin=-50,vmax=50), extent=[ra_relative[0],ra_relative[-1],dec_relative[0],dec_relative[-1]])
ax1 = plt.gca()
ax1.set_aspect('equal')

plt.xlabel(r'RA - RA$\mathrm{_{centre}}$ (arcsec)',fontsize=15)
plt.ylabel(r'Dec. - Dec.$\mathrm{_{centre}}$ (arcsec)',fontsize=15)
ax1.tick_params(axis='y', which='both', right=False, ); ax1.tick_params(axis='x', which='both', top=False )
secax2 = ax1.secondary_xaxis('top', functions=(ang2pctrans, ang2pctrans_inv))
secax2.set_xticklabels([])
secax = ax1.secondary_yaxis('right', functions=(ang2pctrans, ang2pctrans_inv))
secax.set_ylabel(r'Dec. - Dec.$\mathrm{_{centre}}$ (pc)',fontsize=15,rotation=270,labelpad=20) # edited by Eric

ax1.contour(ra_relative_cloud, dec_relative_cloud, smooth_cloud, levels=[7.5e-3,3.0e-2, 6.5e-2, 0.1], colors='purple', zorder=1)
# ,alpha=0.4

scalebar(ax1,)

plt.text(11,10.4,'$HST$',fontsize=11)
plt.text(11,9.2,'ALMA CO(2-1)',color='purple',fontsize=11)

plt.savefig('/Users/ericliang/n1387/work_pub/plot/hst_lower.pdf') # and hst_lower-original_kernel.pdf and other colours
plt.show()




''' molecular mass calculation '''

lines = open('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat').readlines()

flux_total = 83.3 # Jy km/s
error = 8.3
alpha_co = 4.3 # alpha_10, Bolatto+13
frequency_intrinsic = 115.2712 # GHz
distance = float(lines[14].strip()) / 1e6 # MPc
z_helio = float(lines[24].strip())
inc = float(lines[6].strip())

print(distance); print(z_helio)

mass = 3.25e7 * alpha_co / frequency_intrinsic**2 / (1+z_helio)  * flux_total * distance**2
mass = 1.05e4 * (alpha_co/4.3) * flux_total * distance**2 / (1+z_helio) # Bolatto eq 3
print('Molecular gas mass %.3e M_sun' %(mass))
print('error %.1e M_sun' %(error/flux_total*mass))
print('Molecular gas mass in log %.4f M_sun' %(np.log10(mass)))
print('error %.4f M_sun' %( np.log10( mass + error/flux_total*mass) - np.log10(mass) ))
print('Fraction of stellar mass:', mass/4.7e10)

radius_disc = 900 # estimated by eye, pc
surface_density = mass / (np.pi * radius_disc**2) # no factor of np.cos(np.deg2rad(inc))
print(surface_density) # 104 M_sun / pc^2

''' ALMA spectral coverage '''

lines = open('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat').readlines()
f0 = 230.538
z = float(lines[24].strip())
f_center = f0 / (1+z)
width = 977e-6 # 2.0 # 488e-6 # 1.875 # GHz
f_up = f_center + width/2; f_low = f_center - width/2
f_up_gal = f_up * (1+z); f_low_gal = f_low * (1+z)
v_up_gal = 299792.458 * (1-f_up_gal/f0); v_low_gal = 299792.458 * (1-f_low_gal/f0) # radio convention
v_width = v_low_gal-v_up_gal
print('%.4f km/s' %v_width)


''' moment maps '''

import pymakeplots
from astropy.coordinates import Angle

dir_cube = '/Users/ericliang/n1387/work_pub/data/'
cube_corr = dir_cube+'NGC1387_combine_clean5.image.pbcor.fits'
cube_flat = dir_cube+'NGC1387_combine_clean5.image.fits'
# mask = dir_cube+'NGC1387_combine_clean9.mask.fits'
lines = open('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat').readlines()

''' only using RMS mask 202405 '''

plotter=pymakeplots.pymakeplots(cube=cube_corr,cube_flat=cube_flat, clean_mask=None)

ra = Angle(lines[10].strip()); dec = Angle(lines[12].strip()); plotter.obj_ra = ra.degree; plotter.obj_dec = dec.degree
plotter.posang = float(lines[4].strip()); plotter.inc=float(lines[6].strip())
plotter.gal_distance = float(lines[14].strip()) / 1e6; plotter.vsys = float(lines[16].strip()); plotter.redshift = float(lines[24].strip())

plotter.imagesize = [9.999,9.999] # [9.999,9.999] # [10.0,10.0] # half side length in arcsec
plotter.chans2do=[8,93] # emission exists in [a,b), [8,93], testing [2,98]
plotter.mom0_tick = [0,0.03,0.06,0.09,0.12,0.15]; 
plotter.vticks = [-60, -30, 0, 30, 60]; plotter.max_velocity = 80
plotter.maxvdisp = 12 # 12 should be proper

plotter.spatial_smooth=2
plotter.spectral_smooth=2
plotter.rmsfac = 5. # rms clipping
plotter.holesize = -1; plotter.islandsize = -1 # pixels

# store=plotter.make_moments(mom=[3,0,1,2],pdf=True,fits=False,circle=[300,600],error=False,title_text=r'mask$_{\rm RMS}$')
store=plotter.make_moments(mom=[0,1,2],pdf=True,fits=False,circle=[300,600],error=False)


# the below needs input of the flux scaling factor from total_flux.py
# add in 0 for empty channels in NGC1387_mask.fits
scaling = np.array([ 0., 0.46412891, -0.        , -0.        ,  0.21486593,  0.39502385,
        0.52544579,  0.66535705,  0.69005517,  0.6838764 ,  0.71204194,
        0.62932675,  0.62630663,  0.59354837,  0.61546632,  0.65474499,
        0.66466261,  0.69686414,  0.76175659,  0.70806061,  0.73813581,
        0.7145078 ,  0.68664007,  0.79778054,  0.67432729,  0.70200825,
        0.71912262,  0.6315791 ,  0.67810554,  0.73800646,  0.72167869,
        0.70402846,  0.74139565,  0.66589365,  0.65561798,  0.72769008,
        0.66985256,  0.83066873,  0.73647024,  0.67355395,  0.68651213,
        0.62961311,  0.65159637,  0.64155698,  0.67547307,  0.70764152,
        0.70474845,  0.72478209,  0.77228009,  0.66374448,  0.66308304,
        0.6428236 ,  0.69612196,  0.68119225,  0.74166795,  0.75696987,
        0.71330903,  0.7015165 ,  0.68750884,  0.7020529 ,  0.73949594,
        0.74821888,  0.7040332 ,  0.69121984,  0.69824034,  0.67886643,
        0.67534195,  0.67320165,  0.70601604,  0.74384025,  0.67738899,
        0.66908105,  0.62285314,  0.68985357,  0.72833384,  0.72959658,
        0.72540244,  0.71252223,  0.64573862,  0.45415806,  0.30030202,
        0.20996302,  0.10408077, -0.        ,  0.        ])

plotter.range_spec = 100 # showing v_sys pm this value
_ = plotter.make_spec(pdf=True,fits=False,update_v=False,scaling=scaling, dialation=40) 
# dialation=40 with mask, 254.502 +- 2.541 Jy km/s
# v_mean (mask applied): 1286.416 +- 0.442 km/s

# line ratio
r_21 =  255 / 83.3 / (230.538 / 115.2712)**2
r21_err = np.sqrt((3/255)**2 + (8.3/83.3)**2 +0.1**2) * r_21
print('%.4f +- %.2f' %(r_21,r21_err) )  # 0.77 +- 0.08
print('alpha:', 4.3 / r_21) # 5.62


''' continuum map for appendix '''

import pymakeplots
from astropy.coordinates import Angle

dir_cube = '/Users/ericliang/n1387/work_pub/data/'
lines = open(dir_cube+'/NGC1387_galpar.dat').readlines()
cube_corr = dir_cube+'NGC1387_combine_clean5.image.pbcor.fits'
cube_flat = dir_cube+'NGC1387_combine_clean5.image.fits'

plotter=pymakeplots.pymakeplots(cube=cube_corr,cube_flat=cube_flat) #  , clean_mask=None

ra = Angle(lines[10].strip()); dec = Angle(lines[12].strip()); plotter.obj_ra = ra.degree; plotter.obj_dec = dec.degree
plotter.gal_distance = float(lines[14].strip()) / 1e6; 
plotter.imagesize = [9.999,9.999] # [9.999,9.999] # [10.0,10.0] # half side length in arcsec
plotter.chans2do=[8,93] # emission exists in [a,b), [8,93], testing [2,98]
plotter.posang = float(lines[4].strip()); plotter.inc=float(lines[6].strip())

plotter.make_continuum(cont_file=dir_cube+'/NGC1387_continuum_2.image.pbcor.fits',radius=13,vmax=0.4)

corr_slice = fits.getdata(cube_corr)[50]
flat_slice = fits.getdata(cube_flat)[50]
pb = flat_slice / corr_slice
cont_corr = fits.getdata(dir_cube+'/NGC1387_continuum_2.image.pbcor.fits')
cont_hdr = fits.getheader(dir_cube+'/NGC1387_continuum_2.image.pbcor.fits')
cont_flat = pb * cont_corr
fits.writeto(dir_cube+'/NGC1387_continuum_2.image.manual_flat.fits', cont_flat, cont_hdr)

# significance measurement of second northwestern faint diffuse source
# beam 0.149498" X 0.129769", beam size = 1.1331 * bmaj * bmin = 0.022 arcsec^2
# aperture size pi*a*b = np.pi*0.448494*0.389306 = 0.549 arcsec^2 = 24.9 * beam_size
# aperture flux (PB uncorrected) 1.4 mJy
# RMS 0.02 mJy/beam (from CARTA)
# total noise within aperture: 0.02 * np.sqrt(24.9) = 0.10 mJy, 
# S/N=1.4 / 0.10 = 14
# aperture flux (PB corrected) 2.2 mJy
# noise = 2.2/14 = 0.2 mJy

''' continuum map for appendix end '''


''' rotation curve & circular velocity curve '''

plotter.range_spec = 100 # km/s
plotter.pvd_radius = 10 # arcsec
plotter.pvdthick = 1.5 # half-width, unit of arcsec, default=1
plotter.make_pvd(pdf=True,fits=True,vel_file='/Users/ericliang/n1387/work_pub/data/NGC1387_velocity_curve.dat')


''' old draft 202404 '''

plotter=pymakeplots.pymakeplots(cube=cube_corr,cube_flat=cube_flat, clean_mask=mask)
# Factor to convert Jy/beam to K: 1014.981

ra = Angle(lines[10].strip())
dec = Angle(lines[12].strip())
plotter.obj_ra = ra.degree
plotter.obj_dec = dec.degree
plotter.posang = float(lines[4].strip()) # from Hope
plotter.inc=float(lines[6].strip())
plotter.gal_distance = float(lines[14].strip()) / 1e6
plotter.vsys = float(lines[16].strip())
plotter.redshift = float(lines[24].strip())

plotter.imagesize = [10.,10.] # half side length in arcsec # [9.99,9.99] for the images # [10.,10.] for the FITS files
plotter.chans2do=[8,93] # safely away from emission
plotter.rmsfac = 1 # rms clipping
# plotter.maxmom0 = None # 0.15
# plotter.holesize = 100 # now using default 2 beam size
# plotter.islandsize = 16 # now using default 1 beam size

# plotter.make_all(pdf=True,fits=False)

plotter.mom0_tick = [0,0.03,0.06,0.09,0.12,0.15]
plotter.subscript = ''
plotter.max_velocity = 80
plotter.vticks = [-60, -30, 0, 30, 60]
plotter.maxvdisp = 12.

# only set fits and error to True in the final run
store = plotter.make_moments(mom=[0,1,2],pdf=False,fits=True,circle=[300,600],error=True)

plotter.range_spec = 100 # v_sys pm this value
plotter.make_spec(pdf=True,fits=False,update_v=False)
# with mask, 304.157 +- 1.029 Jy km/s
# v_mean (mask applied): 1285.261 km/s
r_21 =  304.157 / 83.3 / (230.538 / 115.2712)**2
r21_err = np.sqrt((1.029/304.157)**2 + (8.3/83.3)**2) * r_21
print('%.4f +- %.2f' %(r_21,r21_err) )

# plotter.make_pvd(pdf=True,fits=True)



''' to show masking process '''
''' Appendix A '''
import pymakeplots
from astropy.coordinates import Angle

dir_cube = '/Users/ericliang/n1387/work_pub/data/'
cube_corr = dir_cube+'NGC1387_combine_clean5.image.pbcor.fits'
cube_flat = dir_cube+'NGC1387_combine_clean5.image.fits'
mask = dir_cube+'NGC1387_combine_clean9.mask.fits'
mask_final = 'NGC1387_mask.fits'

lines = open('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat').readlines()

''' mask_final '''

plotter=pymakeplots.pymakeplots(cube=cube_corr,cube_flat=cube_flat, clean_mask=mask)

ra = Angle(lines[10].strip())
dec = Angle(lines[12].strip())
plotter.obj_ra = ra.degree
plotter.obj_dec = dec.degree
plotter.posang = float(lines[4].strip()) # from Hope
plotter.inc=float(lines[6].strip())
plotter.gal_distance = float(lines[14].strip()) / 1e6
plotter.vsys = float(lines[16].strip())
plotter.redshift = float(lines[24].strip())
plotter.imagesize = [9.99,9.99] # half side length in arcsec
plotter.chans2do=[8,93] # safely away from emission
plotter.rmsfac = 1 # rms clipping
plotter.mom0_tick = [0,0.03,0.06,0.09,0.12,0.15]
plotter.subscript = ''
plotter.max_velocity = 80
plotter.vticks = [-60, -30, 0, 30, 60]
plotter.maxvdisp = 12.


plotter.make_moments(mom=[3,0,1,2],pdf=True,fits=False,circle=[300,600],error=False, title_text=r'mask$_{\rm final}$')
# clean mask sum 1,739,424 --> 863,247
# With mask, total flux in Jy km/s 304.15723


''' clean mask '''

plotter=pymakeplots.pymakeplots(cube=cube_corr,cube_flat=cube_flat, clean_mask=mask)

ra = Angle(lines[10].strip())
dec = Angle(lines[12].strip())
plotter.obj_ra = ra.degree
plotter.obj_dec = dec.degree
plotter.posang = float(lines[4].strip()) # from Hope
plotter.inc=float(lines[6].strip())
plotter.gal_distance = float(lines[14].strip()) / 1e6
plotter.vsys = float(lines[16].strip())
plotter.redshift = float(lines[24].strip())
plotter.imagesize = [9.99,9.99] # half side length in arcsec
plotter.chans2do=[8,93] # safely away from emission
plotter.mom0_tick = [0,0.03,0.06,0.09,0.12,0.15]
plotter.subscript = ''
plotter.max_velocity = 80
plotter.vticks = [-60, -30, 0, 30, 60]
plotter.maxvdisp = 12. # 12.; 1005. for testing

plotter.rmsfac = -1e8 # rms clipping
plotter.holesize = -1 # pixels
plotter.islandsize = -1 # pixels

store = plotter.make_moments(mom=[3,0,1,2],pdf=True,fits=False,circle=[300,600],error=False, title_text=r'mask$_{\rm clean}$')
# clean mask sum 1,739,424 --> 1,739,424
# With mask, total flux in Jy km/s 327.17514

# note Clean mask RA max (arcsec) 10.108345354653107; Clean mask Dec max (arcsec) 10.15909718151704. So, using whole image:
# with mask, total flux in Jy km/s 327.17432 # A difference is expected; being tiny means the process is robust 


''' RMS mask '''

plotter=pymakeplots.pymakeplots(cube=cube_corr,cube_flat=cube_flat, clean_mask=None)

ra = Angle(lines[10].strip())
dec = Angle(lines[12].strip())
plotter.obj_ra = ra.degree
plotter.obj_dec = dec.degree
plotter.posang = float(lines[4].strip()) # from Hope
plotter.inc=float(lines[6].strip())
plotter.gal_distance = float(lines[14].strip()) / 1e6
plotter.vsys = float(lines[16].strip())
plotter.redshift = float(lines[24].strip())
plotter.imagesize = [9.99,9.99] # half side length in arcsec
plotter.chans2do=[8,93] # safely away from emission
plotter.mom0_tick = [0,0.03,0.06,0.09,0.12,0.15]
plotter.subscript = ''
plotter.max_velocity = 80
plotter.vticks = [-60, -30, 0, 30, 60]
plotter.maxvdisp = 30. # special for mask_RMS # using 80 for paper

plotter.spectral_smooth=1
plotter.rmsfac = 1 # rms clipping
plotter.holesize = -1 # pixels
plotter.islandsize = -1 # pixels

store = plotter.make_moments(mom=[3,0,1,2],pdf=True,fits=False,circle=[300,600],error=False, title_text=r'mask$_{\rm RMS}$')
# final mask sum 2,651,610
# With mask, total flux in Jy km/s 802.368. Too much noise picked up above a positive threshold.

# Testing with RMS clipping at 3 sigma:
# With mask, total flux in Jy km/s 112.06078. Missing flux.
# by the way, entire 10"-radius box sum * dv: Without mask, raw total sum flux in plotted region Jy km/s 259.8327 


''' third row, clean * RMS '''
''' change fontsize of mask plot title from 50 to 35 '''

plotter=pymakeplots.pymakeplots(cube=cube_corr,cube_flat=cube_flat, clean_mask=mask)

ra = Angle(lines[10].strip())
dec = Angle(lines[12].strip())
plotter.obj_ra = ra.degree
plotter.obj_dec = dec.degree
plotter.posang = float(lines[4].strip()) # from Hope
plotter.inc=float(lines[6].strip())
plotter.gal_distance = float(lines[14].strip()) / 1e6
plotter.vsys = float(lines[16].strip())
plotter.redshift = float(lines[24].strip())
plotter.imagesize = [9.99,9.99] # half side length in arcsec
plotter.chans2do=[8,93] # safely away from emission
plotter.mom0_tick = [0,0.03,0.06,0.09,0.12,0.15]
plotter.subscript = ''
plotter.max_velocity = 80
plotter.vticks = [-60, -30, 0, 30, 60]
plotter.maxvdisp = 12. # special for mask_RMS

plotter.rmsfac = 1 # rms clipping
plotter.holesize = -1 # pixels
plotter.islandsize = -1 # pixels

store = plotter.make_moments(mom=[3,0,1,2],pdf=True,fits=False,circle=[300,600],error=False, title_text=r'mask$_{\rm{clean}}$$~\times$ mask$_{\rm{RMS}}$')
# clean mask sum 1,739,424 --> 861,820
# With mask, total flux in Jy km/s 305.62686







# plotter=pymakeplots.pymakeplots(cube=cube_corr,cube_flat=cube_flat) # for "mask_rms" version
# # ... #
# plotter.maxvdisp = 80. # for "mask_rms" version

# plotter.rmsfac = -100000 # cancel rms clipping, for "mask_clean" version

# plotter.dpi = 100
# store = plotter.make_moments(mom=[0,1,2],pdf=True,fits=False,circle=[300,600],error=False)


''' old plotting of all masks '''
# mask_clean = fits.getdata('/Users/ericliang/n1387/work_pub/output_publication/maps-v15-appendix/clean/NGC1387_mask_mom.fits')
# # mask_clean = mask_clean[:,plotter.spatial_trim[0]:plotter.spatial_trim[1],plotter.spatial_trim[2]:plotter.spatial_trim[3]]

# mask_rms = fits.getdata('/Users/liangf/work/output_publication/maps-v15-appendix/RMS/NGC1387_mask_mom.fits')

# mask_inter = fits.getdata('/Users/liangf/work/output_publication/maps-v15-appendix/intermediate/NGC1387_mask_mom.fits')

# mask_final = fits.getdata('/Users/liangf/work/output_publication/maps-v15-appendix/final/NGC1387_mask_mom.fits')


# fig,ax = plt.subplots(4,1,figsize=(3,14))
# plt.sca(ax[0])
# plt.imshow(mask_clean[17])
# plt.title('Mask$_{\\rm clean}$')
# plt.sca(ax[1])
# plt.imshow(mask_rms[17])
# plt.title('Mask$_{\\rm RMS}$')
# plt.sca(ax[2])
# plt.imshow(mask_inter[17])
# plt.title('Mask$_{\\rm clean}~\\times$ Mask$_{\\rm RMS}$')
# plt.sca(ax[3])
# plt.imshow(mask_final[17])
# plt.title('Combined mask after \n removal of holes and islands')
# ax[0].axis('off');ax[1].axis('off');ax[2].axis('off');ax[3].axis('off')
# # plt.tight_layout()
# # plt.subplots_adjust()
# # plt.show()
# plt.savefig('masking.pdf') # ,bbox_inches='tight'
''' old plotting of all masks END'''


''' two remaining possibilities not shown: clean+removal, rms+removal '''




''' testing orientation of beam '''
# plotter.bmaj = plotter.bmaj * 50
# plotter.bmin = plotter.bmin / 10 
# plotter.bpa = 70
# store = plotter.make_moments(mom=[0,],pdf=True,fits=True,circle=[300,600],error=False)

''' flux difference with Pandora's cube '''
import pymakeplots
from astropy.coordinates import Angle

dir_cube = '/Users/ericliang/n1387/work_pub/data/pandora/'
cube_corr = dir_cube+'NGC1387.CO.0.8-high_res-pandora.fits'

plotter_p =pymakeplots.pymakeplots(cube=cube_corr)
# Factor to convert Jy/beam to K: 1734.193

beam_ratio = (0.167104*0.135577)/(0.128863*0.102899) # 1.7, for CARTA viewing

lines = open('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat').readlines()
ra = Angle(lines[10].strip())
dec = Angle(lines[12].strip())
plotter_p.obj_ra = ra.degree
plotter_p.obj_dec = dec.degree
plotter_p.posang = float(lines[4].strip()) # from Hope
plotter_p.inc=float(lines[6].strip())
plotter_p.gal_distance = float(lines[14].strip()) / 1e6
plotter_p.vsys = float(lines[16].strip())
plotter_p.redshift = float(lines[24].strip())

plotter_p.imagesize = [9.99,9.99] # half side length in arcsec
plotter_p.chans2do=[8,93] # safely away from emission
plotter_p.rmsfac = 1 # rms clipping

...


''' CPROPS flux fraction '''

cube = fits.getdata('/Users/liangf/work/data/NGC1387_combine_clean5.image.pbcor.fits')
mask_c = fits.getdata('/Users/liangf/work/measurements_publication/NGC1387_CO21_cube_2kms_mask.fits')
mask_clean = fits.getdata('/Users/liangf/work/output_publication/maps-v15-appendix/clean/NGC1387_mask_mom.fits')
mask_moment = fits.getdata('/Users/liangf/work/output_publication/maps-v15-appendix/final/NGC1387_mask_mom.fits')

dv = 2
beam_per_spaxel = 1/16.
flux_c = np.nansum(cube*mask_c) * dv * beam_per_spaxel

spatial_trim = plotter.spatial_trim
chans2do = plotter.chans2do
mask_c_trim = mask_c[chans2do[0]:chans2do[1], spatial_trim[2]:spatial_trim[3], spatial_trim[0]:spatial_trim[1]]
cube_trim = cube [chans2do[0]:chans2do[1], spatial_trim[2]:spatial_trim[3], spatial_trim[0]:spatial_trim[1]]

flux_clean = np.sum(cube_trim * mask_clean) * dv * beam_per_spaxel
flux_moment = np.sum(cube_trim * mask_moment) * dv * beam_per_spaxel
flux_raw = np.sum(cube_trim) * dv * beam_per_spaxel

# these are missing dv=2 and beam_area / spaxel_area
print('Raw sum: %.2f Jy' %flux_raw) # 264.07
print('Clean mask sum: %.2f Jy' %flux_clean)  # 328.08
print('Moment mask sum: %.2f Jy' %flux_moment) # 288.96, used in the paper
print('CPROPS mask sum: %.2f Jy' %flux_c) # 259.57




''' Plot the continuum image old '''

continuum_hdu = fits.open('/Users/liangf/work/ngc1387/calibrated_all/continuum_trial/NGC1387_continuum_2.image.pbcor.fits')
continuum = fits.getdata('/Users/liangf/work/ngc1387/calibrated_all/continuum_trial/NGC1387_continuum_2.image.pbcor.fits')

plt.imshow(continuum)
plt.savefig('n1387-continuum.pdf')
plt.show()


# # ACA map & FoV

# filename = '/Users/fuhengliang/ox_local/alma_archive/ngc1387/calibrated_all/NGC_1387_CO_cube.clean2.moment.fits'

# hdu_aca = fits.open(filename)
# hdu = fits.open(filename)
# map_ = hdu[0].data

# center = len(map_)/2
# unit_size = np.absolute(hdu[0].header['CDELT2']) * 3600.
# fov_aca = 45.7 # 39.
# fov_12m = 27.4 # 23.

# fig,ax = plt.subplots(1,1)
# plt.imshow(map_, vmin=1)
# cb = plt.colorbar()
# cb.set_label('Jy/beam * km/s')
# plt.title('NGC1387_mom0 ACA')
# circle1 = plt.Circle((center, center), fov_12m/2./unit_size, facecolor='none', edgecolor='green', ls='--', label='12-m FoV')
# circle2 = plt.Circle((center, center), fov_aca/2./unit_size, facecolor='none', edgecolor='red', ls='--', label='ACA FoV')
# ax.add_patch(circle1)
# ax.add_patch(circle2)
# plt.legend()
# plt.show()



'''
# test the difference of WCS projection

hst_shape = hst_image.shape
origin = SkyCoord.from_pixel(0,0,hst_wcs)
end = SkyCoord.from_pixel(hst_shape[0],hst_shape[1],hst_wcs)
# print('RA increment (all correction): %.5f' %( np.absolute(origin.ra.arcsecond - end.ra.arcsecond)/float(hst_shape[0]) ))
# origin_core = SkyCoord.from_pixel(0,0,hst_wcs,mode='wcs')
# end_core = SkyCoord.from_pixel(hst_shape[0],hst_shape[1],hst_wcs,mode='wcs')
# print('RA increment (core correction): %.5f' %( np.absolute(origin_core.ra.arcsecond - end_core.ra.arcsecond)/float(hst_shape[0]) ))
# No difference, suspect that no extra correction term is provided other than core terms in the header. Therefore, "all" mode cannot do anything more.
# print('DEC increment: %.5f' %( np.absolute(origin.dec.arcsecond - end.dec.arcsecond)/float(hst_shape[1]) ))

# end test

ra_increment = np.absolute(origin.ra.arcsecond - end.ra.arcsecond)/float(hst_shape[0])
dec_increment = np.absolute(origin.dec.arcsecond - end.dec.arcsecond)/float(hst_shape[0])

cloud_hdu = fits.open('/Users/liangf/work/ngc1387/calibrated_all/NGC1387_mom0.fits')
cloud_map = cloud_hdu[0].data
kernel = Gaussian2DKernel(x_stddev=4)
smooth_cloud = convolve(cloud_map, kernel)

gal_distance=18.5; barlength_pc = 1000; barlength_arcsec=barlength_pc / (gal_distance*1e6) / np.pi * 180. * 3600.
pix_number = barlength_arcsec / ra_increment
gal_distance_upp = 19.5; error=pix_number*ra_increment * (gal_distance_upp*1e6) * np.pi / 180. / 3600. / 1000 - 1
label=r'1$\pm$%.2f kpc'%error

center = SkyCoord(54.237750, -35.506639,unit='deg',frame='fk5')
center_pixel = center.to_pixel(hst_wcs)
ra_side = 220; dec_side = ra_side * ra_increment / dec_increment

x_low = int(center_pixel[0]-ra_side); x_high = int(center_pixel[0]+ra_side)
y_low = int(center_pixel[1]-dec_side); y_high = int(center_pixel[1]+dec_side)

# hst_image_norm = hst_image / np.max(hst_image)
# enhanced = unsharp_mask(hst_image_norm, radius=20, amount=2) * np.max(hst_image)
# test
enhanced = unsharp_mask(hst_image, radius=5, amount=1, preserve_range=True)
difference  = enhanced - hst_image
blurred = -1 * (enhanced - 2 * hst_image)

hst_median = np.zeros(hst_image.shape)
hst_median[ y_low : y_high , x_low : x_high ] = ndimage.median_filter(hst_image[y_low:y_high,x_low:x_high], size=20)
difference_median = hst_image - hst_median
dust = difference_median[ y_low : y_high , x_low : x_high ] ** 2 # external

plt.figure()
ax = plt.subplot(projection=hst_wcs)
ax_ra = ax.coords[0]; ax_dex = ax.coords[1]
ax_ra.set_major_formatter('dd:mm:ss.ss')

# plt.imshow(hst_image, origin='lower', aspect=dec_increment/ra_increment, cmap='inferno', norm=LogNorm(vmin=0.2,vmax=50)) # #,extent=[1,2,3,4] #extent=[54.24167, 54.233333, -35.51083, -35.5025]
# plt.imshow(difference, origin='lower', aspect=dec_increment/ra_increment, cmap='inferno', vmin=-3, vmax=3) # #,extent=[1,2,3,4] #extent=[54.24167, 54.233333, -35.51083, -35.5025] , norm=LogNorm(vmin=0.2,vmax=70) cm.vik norm=SymLogNorm(linthresh=0.5, vmin=-3,vmax=3)
# plt.title('Enhanced = Original - Blurred')

plt.imshow(difference_median, origin='lower', aspect=dec_increment/ra_increment, cmap='inferno', vmin=-3, vmax=3) # #,extent=[1,2,3,4] #extent=[54.24167, 54.233333, -35.51083, -35.5025] , norm=LogNorm(vmin=0.2,vmax=70) cm.vik   , norm=SymLogNorm(linthresh=0.5, vmin=-3,vmax=3)
plt.title('Enhanced = Original - Median')

plt.colorbar()

ax.contour(smooth_cloud, transform=ax.get_transform(wcs.WCS(cloud_hdu[0].header)),
            levels=[7.5e-3,3.0e-2, 6.5e-2, 0.1],colors='white',alpha=0.4) # levels=[1.75e-2,4.58e-2,7.74e-2,0.11]

# ell = Ellipse(xy=(54.2341,-35.50861), height=0.167104/3600., width=0.135577/3600., 
    # angle=89.5607, color='green',transform=ax.get_transform('world'))
# ax.add_patch(ell)

# ell2 = Ellipse(xy=(54.2361,-35.50811), height=0.5/3600., width=0.5/3600., 
#     angle=0, color='green',transform=ax.get_transform('world'))
# ax.add_patch(ell2)
# ell3 = Ellipse(xy=(2900,2450), height=10, width=10, 
#     angle=0, color='white',transform=ax.get_transform('pixel'))
# ax.add_patch(ell3)

fov_aca = 45.7 # 39.
fov_12m = 27.4 # 23.
circle1 = plt.Circle((center.ra.degree, center.dec.degree), fov_12m/2./3600., facecolor='none', edgecolor='green', ls='--', label='12-m FoV %.2f\"'%fov_12m, transform=ax.get_transform('world'))
circle2 = plt.Circle((center.ra.degree, center.dec.degree), fov_aca/2./3600., facecolor='none', edgecolor='red', ls='--', label='ACA FoV %.2f\"'%fov_aca, transform=ax.get_transform('world'))
# ax.add_patch(circle1)
# ax.add_patch(circle2)
plt.legend()

asb = AnchoredSizeBar(transform=ax.transData, size=pix_number, label=label, loc='lower center', pad=0.5, borderpad=0.3, sep=5,frameon=False,color='white')
ax.add_artist(asb)
plt.xlim(x_low,x_high)
plt.ylim(y_low,y_high)

# ax.coords.grid(True, color='white', ls='solid')

plt.xlabel('RA')
plt.ylabel('Dec')

ax2 = ax.twiny()
lim = np.array(ax.get_xlim())
# lim_convert = 18.5 * lim / 3600. / 180. * np.pi * 1e6
# plt.xlim(lim_convert)
plt.xlabel('RA offset [pc]')

plt.show();
plt.savefig('HST_ALMA.pdf');plt.close()
'''