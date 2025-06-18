# coding: utf-8
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import ICRS
from astropy.table import Table
from astropy import wcs # new
from scipy import ndimage
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Ellipse
from matplotlib import cm
from matplotlib.colors import ListedColormap #, LinearSegmentedColormap
from matplotlib.offsetbox import AnchoredText
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse,AnchoredSizeBar
from pafit.fit_kinematic_pa import fit_kinematic_pa
from spectral_cube import SpectralCube
from spectral_cube.utils import SpectralCubeWarning
from skimage import morphology # added by Eric
from uncertainties import unumpy # added by Eric
import warnings
from pymakeplots.sauron_colormap import sauron
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.ndimage import binary_dilation as dia

# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.colors import LogNorm

warnings.filterwarnings(action='ignore', category=SpectralCubeWarning, append=True)
warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning, append=True)
warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)

def running_mean(x, N): # used in make_spec
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def rotateImage(img, angle,): # used in make_PVD pivot
    # padX = [img.shape[0] - pivot[0], pivot[0]]
    # padY = [img.shape[1] - pivot[1], pivot[1]]
    # padZ = [0,0]
    # imgP = np.pad(img, [padY, padX, padZ], 'constant')
    imgR = ndimage.rotate(img, angle, reshape=True, mode='constant', cval=np.nan)
    return imgR


class pymakeplots:
    def __init__(self,cube_flat=None,pb=None,cube=None,clean_mask=None): # added clean_mask by Eric
        self.galname=None
        self.gal_distance=None # Mpc
        self.posang=None
        self.inc=0.
        self.vsys=None
        self.moment1=None
        self.rms=None
        self.flat_cube=None
        self.pbcorr_cube=None
        self.spectralcube=None
        self.mask=None
        self.flat_cube_trim=None
        self.pbcorr_cube_trim=None
        self.mask_trim=None
        self.bmaj=None
        self.bmin=None
        self.bpa=None
        self.xcoord,self.ycoord,self.vcoord = None, None, None
        self.xcoord_trim,self.ycoord_trim,self.vcoord_trim = None, None, None
        self.dv=None
        self.cellsize=None
        self.silent=False # rig for silent running if true
        self.bardist=None
        self.rmsfac=1.5
        self.holesize=None # added by Eric
        self.islandsize=None # added by Eric
        self.spatial_smooth=1.0 # in unit of beam FWHM (major axis)
        self.spectral_smooth=1.0 # in unit of channels
        self.restfreq=None
        self.obj_ra=None
        self.obj_dec=None
        self.imagesize=None
        self.xc=None
        self.yc=None
        self.bunit=None
        # self.linefree_chans_start, self.linefree_chans_end = 0, 6 # edited by Eric
        self.chans2do=None
        self.spatial_trim = None
        self.maxvdisp=None
        self.cliplevel=None
        self.fits=False
        self.pvdthick=1.  # half-width, unit of arcsec, modified by Eric
        self.flipped=False # for possible mismatch of velocity/frequency axis between different inputs
        self.make_square=True
        self.useallpixels = False
        self.wcs=None
        self.clean_mask=None
        self.pb = None # added by Eric
        self.pb_trim = None # added by Eric
        self.dpi = None # added by Eric
        self.mom0_tick = None # added by Eric
        self.maxmom0 = None # added by Eric
        self.plot_axis = True # added by Eric Kinematic major axis
        self.range_spec= None # added by Eric
        self.pvd_radius= None
        self.vticks = None
        self.subscript=''
        self.max_velocity = None # for moment 1 map
        self.centre_marker = False
        self.redshift = 0.0 
        self.cb_pad = 0.09
        # self.cb_frac = 0.15 # default of matplotlib
        
        # every cube read in below is transposed by read_in_a_cube()

        if (cube != None)&(pb==None)&(cube_flat==None):
            # only one cube given
            self.input_cube_nopb(cube)
            
        if (cube != None)&(pb==None)&(cube_flat!=None):
            # pbcorred cube and flat cube given
            self.input_cube_pbcorr_and_flat(cube,cube_flat)

        if (cube != None)&(pb!=None):
            # pbcorred cube and pb given, whether or not cube_flat is given
            self.input_cube_pbcorr(cube,pb)
        
        if (cube ==None)&(pb!=None)&(cube_flat!= None): # edited by Eric
            # flat cube and pb given
            self.input_cube_flat(cube_flat,pb)  
              
        # check, added by Eric
        if np.any(self.pbcorr_cube) == None:
            print('Error: missing essential input')

        # added by Eric
        self.rms = self.rms_estimate(self.flat_cube) 
        print('\n RMS in unsmoothed uncorrected cube (centre of channel #1) %.4e %s' %(self.rms, self.bunit))
        self.pb = self.flat_cube[:,:,0] / self.pbcorr_cube[:,:,0]
            
    def vsystrans_inv(self,val):
        return val +self.vsys

    def vsystrans(self,val):
        return val - self.vsys
        
    def ang2pctrans_inv(self,val):
        return val/(1. / 3600. / 180. * np.pi * 1e6 *self.gal_distance)

    def ang2pctrans(self,val):
        return val*1. / 3600. / 180. * np.pi * 1e6*self.gal_distance
        
    def ang2kpctrans_inv(self,val):
        return val/(1. / 3600. / 180. * np.pi * 1e3*self.gal_distance)

    def ang2kpctrans(self,val):
        return val*1. / 3600. / 180. * np.pi * 1e3*self.gal_distance    

    def beam_area(self):
        return (np.pi*(self.bmaj/self.cellsize)*(self.bmin/self.cellsize))/(4*np.log(2))

    def rms_estimate(self,cube,): # edited by Eric
        quarterx=np.array(cube.shape[0]/3.).astype(int)
        quartery=np.array(cube.shape[1]/3.).astype(int)
        return np.nanstd(cube[quarterx*1:2*quarterx,1*quartery:2*quartery,1])

    def pri_beam(self, amp, mean_x, mean_y, fwhm, size_x, size_y):
        ''' 
        mean_x, mean_y, fwhm in unit of arcsec
        size_x, size_y in unit of pixel
        '''
        x, y = np.meshgrid(np.linspace(0.5-size_x/2, 0.5-size_x/2+size_x-1,size_x), np.linspace(0.5-size_y/2, 0.5-size_y/2+size_y-1,size_y))
        d =  (x-mean_x/self.cellsize)**2 + (y-mean_y/self.cellsize)**2
        sigma = fwhm / self.cellsize / (2*np.sqrt(2*np.log(2)))
        beam_array = amp * np.exp(-( d / ( 2.0 * sigma**2 ) ) )
        return beam_array

    def set_rc_params(self,mult=1):
        #matplotlib.rcParams['text.usetex'] = True
        #matplotlib.rcParams['font.family'] = 'Latin Modern Roman'
        matplotlib.rcParams.update({'font.size': 25*mult})
        matplotlib.rcParams['legend.fontsize'] = 17.5*mult
        matplotlib.rcParams['axes.linewidth'] = 1.5
        matplotlib.rcParams['xtick.labelsize'] = 20*mult
        matplotlib.rcParams['ytick.labelsize'] = 20*mult
        matplotlib.rcParams['xtick.major.size'] = 10
        matplotlib.rcParams['ytick.major.size'] = 10
        matplotlib.rcParams['xtick.major.width'] = 2
        matplotlib.rcParams['ytick.major.width'] = 2
        matplotlib.rcParams['xtick.minor.size'] = 5
        matplotlib.rcParams['ytick.minor.size'] = 5
        matplotlib.rcParams['xtick.minor.width'] = 1
        matplotlib.rcParams['ytick.minor.width'] = 1
        matplotlib.rcParams['xtick.direction'] = 'in'
        matplotlib.rcParams['ytick.direction'] = 'in'
        matplotlib.rcParams['xtick.bottom'] = True
        matplotlib.rcParams['ytick.left'] = True
        matplotlib.rcParams["xtick.minor.visible"] = True
        matplotlib.rcParams["ytick.minor.visible"] = True
        #params = {'mathtext.default': 'regular'}
        #matplotlib.rcParams.update(params)
        matplotlib.rcParams['axes.labelsize'] = 27*mult


        
    def input_cube_nopb(self,path_to_cube):
       
       self.pbcorr_cube = self.read_primary_cube(path_to_cube) 
       
       self.flat_cube = self.pbcorr_cube.copy() # edited by Eric
       
    def input_cube_pbcorr(self,path_to_pbcorr_cube,path_to_pb):
       
       self.pbcorr_cube = self.read_primary_cube(path_to_pbcorr_cube) 
       
       pb,hdr,_= self.read_in_a_cube(path_to_pb)
       if self.flipped: pb=np.flip(pb,axis=2)

       self.flat_cube = self.pbcorr_cube*pb
       self.flat_cube[self.flat_cube == 0.0] = np.nan

       if self.flipped: cube_out = np.flip(self.flat_cube,axis=2)

       hdu = fits.PrimaryHDU(cube_out.T,header = hdr)
       hdul = fits.HDUList([hdu])
       hdul.writeto(self.galname+'_flat.fits', overwrite=True)

    def input_cube_flat(self,path_to_flat_cube,path_to_pb):
       
       self.flat_cube = self.read_primary_cube(path_to_flat_cube)
       
       pb,hdr,_= self.read_in_a_cube(path_to_pb)
       if self.flipped: pb=np.flip(pb,axis=2)
       
       self.pbcorr_cube = self.flat_cube.copy()*0.0
       self.pbcorr_cube[np.isfinite(pb) & (pb != 0)] = self.flat_cube[np.isfinite(pb) & (pb != 0)] / pb[np.isfinite(pb) & (pb != 0)]

    def input_cube_pbcorr_and_flat(self,path_to_pbcorr_cube,path_to_flat_cube):
       
       self.pbcorr_cube = self.read_primary_cube(path_to_pbcorr_cube) 
       
       self.flat_cube,hdr,_ = self.read_in_a_cube(path_to_flat_cube)
       if self.flipped: self.flat_cube=np.flip(self.flat_cube,axis=2)           


    def read_primary_cube(self,cube):
        ### read in cube ###
        datacube,hdr,beamtab = self.read_in_a_cube(cube)
        
        try:
           self.bmaj=np.median(beamtab['BMAJ'])
           self.bmin=np.median(beamtab['BMIN'])
           self.bpa=np.median(beamtab['BPA'])
        except:     
           self.bmaj=hdr['BMAJ']*3600.
           self.bmin=hdr['BMIN']*3600.
           self.bpa=hdr['BPA'] # edited by Eric

        print('BMAJ(") X BMIN(") = %0.4f X %0.4f at %.2f deg' %(self.bmaj, self.bmin, self.bpa))
    
           
        try:
            self.galname=hdr['OBJECT']
        except:
            self.galname="Galaxy"
            
        try:
            self.bunit=hdr['BUNIT']
        except:
            self.bunit="Unknown"
                          
        self.xcoord, self.ycoord, self.vcoord, self.cellsize, self.dv = self.get_header_coord_arrays(hdr,cube_path = cube)

        bmaj = self.bmaj * u.arcsec; bmin = self.bmin * u.arcsec
        beam_area = 1.1330900354567983 * (bmaj*bmin)
        factor = ( (u.Jy/beam_area).to(u.K, equivalencies=u.brightness_temperature(self.restfreq * u.Hz / (1+self.redshift)))).value 
        print('\n Factor to convert Jy/beam to K: %.3f' %factor)
        print('To check, unit in cube:', self.bunit)
        print('To check, unit of frequency being Hz:', hdr.comments['RESTFRQ'])

        # This is run when the instance is created; self.obj_ra and self.obj_dec can still be updated later; trimming is triggered at any plotting task and thus uses the updated ones
        try: 
            self.obj_ra=hdr['OBSRA']
            self.obj_dec=hdr['OBSDEC']
            if (self.obj_ra > np.max(self.xcoord)) or (self.obj_ra < np.min(self.xcoord)) or (self.obj_dec < np.min(self.ycoord)) or (self.obj_dec > np.max(self.ycoord)):
                # obsra/dec given in the headers arent in the observed field! Fall back on medians.
                if not self.silent: 
                    print("OBSRA/OBSDEC keywords dont seem correct! Assuming galaxy centre is at pointing centre")
                    self.obj_ra=np.median(self.xcoord)
                    self.obj_dec=np.median(self.ycoord)
        except:
            self.obj_ra=np.median(self.xcoord)
            self.obj_dec=np.median(self.ycoord)
        
        if self.dv < 0:
            datacube = np.flip(datacube,axis=2)
            self.dv*=(-1)
            self.vcoord = np.flip(self.vcoord)
            self.flipped=True
        datacube[~np.isfinite(datacube)]=0.0
        
        # self.rms= self.rms_estimate(datacube,self.linefree_chans_start,self.linefree_chans_end)  # edited by Eric
        return datacube

    def read_in_a_cube(self,path):
        hdulist=fits.open(path)
        hdr=hdulist[0].header
        cube = np.squeeze(hdulist[0].data.T) #squeeze to remove singular stokes axis if present
        # cube[np.isfinite(cube) == False] = 0.0 # commented out by Eric. Converting nan values to zeros mess up RMS calculation, for example
        
        try:
            if hdr['CASAMBM']:
                beamtab = hdulist[1].data
        except:
            beamtab=None
            
        return cube, hdr, beamtab

    def get_header_coord_arrays(self,hdr,cube_path): # fully rewritten by Eric

        self.wcs=wcs.WCS(hdr)
        print('\n check square pixel, check ra dec order: \n')
        print(self.wcs)

        self.spectralcube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s, velocity_convention='radio') # , 
        try:
            self.restfreq = hdr['RESTFRQ'] # Hz
        except:
            pass

        y,x=self.spectralcube.spatial_coordinate_map

        x1=np.median(x[0:hdr['NAXIS2'],0:hdr['NAXIS1']],0).value
        y1=np.median(y[0:hdr['NAXIS2'],0:hdr['NAXIS1']],1).value
        v1=self.spectralcube.spectral_axis.value

        cd3= np.median(np.diff(v1))
        cd2= np.median(np.diff(y1))
        
        if np.any(np.diff(x1) > 359):
            # we have RA=0 wrap issue
            x1[x1 > 180]-=360
            
            
        return x1,y1,v1,np.abs(cd2*3600),cd3                    


        # zp = np.arange(hdr['NAXIS3'])
        # xp = np.ones(len(zp)) * int(hdr['NAXIS1']/2)
        # yp = np.ones(len(zp)) * int(hdr['NAXIS2']/2)
        # spectral = self.wcs.all_pix2world(xp,yp,zp, 0)[2]

        # xp = np.arange(hdr['NAXIS1'])
        # yp = np.ones(len(xp)) * int(hdr['NAXIS2']/2)
        # zp = np.ones(len(xp)) * int(hdr['NAXIS3']/2)
        # ra = self.wcs.all_pix2world(xp,yp,zp, 0)[0]

        # yp = np.arange(hdr['NAXIS2'])
        # xp = np.ones(len(yp)) * int(hdr['NAXIS1']/2)
        # zp = np.ones(len(yp)) * int(hdr['NAXIS3']/2)
        # dec = self.wcs.all_pix2world(xp,yp,zp, 0)[1]

        # try:
        #     self.restfreq = hdr['RESTFRQ']*u.Hz
        #     print('\n check if unit is Hz \n')
        #     print(hdr.comments['RESTFRQ'])
        # except:
        #     self.restfreq = hdr['RESTFREQ']*u.Hz
        #     print('\n check if unit is Hz \n')
        #     print(hdr.comments['RESTFREQ'])

        # if hdr['CUNIT3'] == 'Hz':
        #     f1=spectral*u.Hz
        #     v1=f1.to(u.km/u.s, equivalencies=u.doppler_radio(self.restfreq))
        #     v1=v1.value
        # elif hdr['CUNIT3'] == 'MHz':
        #     f1=spectral*u.MHz
        #     v1=f1.to(u.km/u.s, equivalencies=u.doppler_radio(self.restfreq))
        #     v1=v1.value
        # elif hdr['CUNIT3'] == 'GHz':
        #     f1=spectral*u.GHz
        #     v1=f1.to(u.km/u.s, equivalencies=u.doppler_radio(self.restfreq))
        #     v1=v1.value
        # elif hdr['CUNIT3']=='m/s':
        #     v1 = spectral/1e3
        # elif hdr['CUNIT3']=='km/s':
        #     v1 = spectral
        # else:
        #     print('No spectral unit found')

        # if np.any(np.diff(ra) > 359):
        #     # we have RA=0 wrap issue
        #     ra[ra > 180]-=360

        # cd3= np.median(np.diff(v1))
        # cd2= np.median(np.diff(dec)) # added by Eric

        # return ra,dec,v1,cd2*3600,cd3 # corrected by Eric


    def prepare_cubes(self):
        
        self.clip_cube()
            
        self.xc=(self.xcoord_trim-self.obj_ra) * 3600 * np.cos(np.deg2rad(self.obj_dec)) # removed by Eric: *(-1)
        self.yc=(self.ycoord_trim-self.obj_dec) * 3600. # removed by Eric: / np.cos(np.deg2rad(self.obj_dec))

        self.mask_trim=self.smooth_mask(self.flat_cube_trim)
        
        if self.gal_distance == None:
            self.gal_distance = self.vsys/70.
            if not self.silent: print("Warning! Estimating galaxy distance using Hubble expansion (h=0.7). Set `gal_distance` if this is not appropriate.")


        if self.posang==None:
            # try fitting the moment one to get the kinematic pa
            if not self.silent: print("No position angle given, estimating using the observed moment one.")

            # mom0=(self.pbcorr_cube_trim*self.mask_trim).sum(axis=2)
            # mom1=mom0.copy()*np.nan
            # mom1[mom0 > 0.0] = (((self.pbcorr_cube_trim*self.mask_trim)*self.vcoord_trim).sum(axis=2))[mom0 > 0.0]/mom0[mom0 > 0.0] # edited by Eric

            cube_temp = self.pbcorr_cube_trim.copy()
            cube_temp[cube_temp<=0] = np.nan # self.cliplevel
            mom0=np.nansum(cube_temp*self.mask_trim, axis=2) # no need for dv here
            mom1=mom0.copy()*np.nan
            mom1 = np.nansum(cube_temp*self.mask_trim*self.vcoord_trim, axis=2) / mom0 # edited by Eric
            mom1[mom1<=0] = np.nan; mom1[~np.isfinite(mom1)] = np.nan
            
            mom1=mom1.T

            # if the cube is small, use it directly to estimate posang. If its large, then interpolate down to keep runtime low.
            if (self.pbcorr_cube_trim[:,:,0].size < 50*50) or (self.useallpixels):
                xv, yv = np.meshgrid(self.xc,self.yc)
                x,y,v = xv[np.isfinite(mom1)],yv[np.isfinite(mom1)],mom1[np.isfinite(mom1)]
                plt.close('all')
                plt.imshow(mom1);plt.show()
            else:
                print("Downsampling the observed moment one in PA estimate for speed. Set `useallpixels` to override.")    
                ind = np.where(np.isfinite(mom1) == False)
                mom1[ind] = self.vsys
                interper = interpolate.interp2d(self.xc,self.yc,mom1-self.vsys,bounds_error=False,fill_value=np.nan)
                # interper = interpolate.LinearNDInterpolator(?,?,mom1-self.vsys,) # bounds_error=False,  fill_value=np.nan
                x=np.linspace(np.min(self.xc),np.max(self.xc),50)
                y=np.linspace(np.min(self.yc),np.max(self.yc),50)
                v= interper(x,y)
                xv, yv = np.meshgrid(x,y)
                x,y,v = xv.flatten(),yv.flatten(),v.flatten()
                
            self.posang, posang_error, velocity_difference = fit_kinematic_pa(-1 * x[np.isfinite(v)],y[np.isfinite(v)],v[np.isfinite(v)],nsteps=36,plot=False,quiet=True)
            print('PA (raw value from fit_kinematic_pa) = %0.2f +- %0.2f'%(self.posang, posang_error))

            mom1 = np.flip(mom1,axis=0)
            # plt.close('all');plt.figure()
            # plt.imshow(mom1,origin='upper');plt.show()
            rotmom1= rotateImage(mom1, -1*self.posang,) # in a frame with "origin=lower", a positive rotation is counter-clockwise
            rotmom1 = np.flip(rotmom1,axis=0)
            # plt.close('all');plt.figure()
            # plt.imshow(rotmom1,origin='upper');plt.show()

            up_median = np.nanmedian(rotmom1[int(rotmom1.shape[0]/2):,:])
            low_median = np.nanmedian(rotmom1[:int(rotmom1.shape[0]/2),:])
            if up_median < low_median:
                self.posang = (self.posang + 180) % 360

            
            # if np.sin(np.deg2rad((self.posang+45)*2)) > 0:
            #     # do y axis cut
            #     if np.nanmean(mom1[self.yc > 0,:]) > np.nanmean(mom1[self.yc < 0,:]):
            #         # posang should be gt 180
            #         if self.posang < 180: self.posang += 180
            #     else:
            #          # posang should be lt 180
            #         if self.posang > 180: self.posang -= 180    
            # else:
            #     # do x axis cut
            #     if np.nanmean(mom1[:,self.xc > 0]) > np.nanmean(mom1[:,self.xc < 0]):
            #         # posang should be gt 180
            #         if self.posang < 180: self.posang += 180
            #     else:
            #          # posang should be lt 180
            #         if self.posang > 180: self.posang -= 180


            # self.posang = 180 - self.posang # edited by Eric, don't fully understand
            # if self.posang < 0:
            #     self.posang += 360
            # if self.posang >= 360:
            #     self.posang -= 360

            if not self.silent: print("PA estimate (degrees): ",np.round(self.posang,1))        


    def clip_cube(self):
        
        # chans2do update begin
        if self.chans2do == None:
            if self.clean_mask is not None: # use the mask to try and guess the channels with signal.
                hdu_clean = fits.open(self.clean_mask)
                clean_mask_data = np.squeeze(hdu_clean[0].data.astype(bool))
                if self.flipped: clean_mask_data=np.flip(clean_mask_data,axis=0)
                index_spec = np.where(clean_mask_data == True)[0] # 0 v, 1 dec, 2 ra
                self.chans2do=[np.clip(np.min(index_spec)-2,0,self.vcoord.size), np.clip(np.max(index_spec)+3,0,self.vcoord.size)]
            else:
                mask_cumsum=np.nancumsum((self.pbcorr_cube > self.rmsfac*self.rms).sum(axis=0).sum(axis=0))
                w_low,=np.where(mask_cumsum/np.max(mask_cumsum) < 0.02)
                w_high,=np.where(mask_cumsum/np.max(mask_cumsum) > 0.98)
                
                if w_low.size ==0: w_low=np.array([0])
                if w_high.size ==0: w_high=np.array([self.vcoord.size])
                self.chans2do=[np.clip(np.max(w_low)-2,0,self.vcoord.size), np.clip(np.min(w_high)+3,0,self.vcoord.size)]
        # chans2do end

        noise_cube_big = self.flat_cube[:,:, np.r_[0:self.chans2do[0],self.chans2do[1]:self.flat_cube.shape[2]]]
        # self.rms = np.nanstd(noise_cube_big)
        self.rms = np.sqrt(np.nanmean(noise_cube_big**2))
        mean_this = np.nanmean(noise_cube_big)
        rms_mean = np.nanstd(noise_cube_big)
        print('In uncorrected cube (all line-free channels), mean = %.4e, RMS (around 0, in use, updated) = %.4e, Std (RMS around mean) = %.4e %s \n' %(mean_this, self.rms, rms_mean, self.bunit))
        
    
        if self.vsys is None:
            # use the cube to try and guess the vsys
            self.vsys=((self.pbcorr_cube*(self.pbcorr_cube > self.rmsfac*self.rms)).sum(axis=0).sum(axis=0)*self.vcoord).sum()/((self.pbcorr_cube*(self.pbcorr_cube > self.rmsfac*self.rms)).sum(axis=0).sum(axis=0)).sum()

        if self.clean_mask is not None:
            hdu_clean = fits.open(self.clean_mask)
            clean_mask_data = np.squeeze(hdu_clean[0].data.astype(bool))
            if self.flipped: clean_mask_data=np.flip(clean_mask_data,axis=0)
            index_dec = np.where(clean_mask_data == True)[1]
            index_ra = np.where(clean_mask_data == True)[2]
            ra_max = np.max(np.absolute( self.xcoord[index_ra] - self.obj_ra)*3600. * np.cos(np.deg2rad(self.obj_dec)))
            dec_max = np.max(np.absolute( self.ycoord[index_dec] - self.obj_dec)*3600.)
            print('Clean mask RA max (arcsec)', ra_max)
            print('Clean mask Dec max (arcsec)', dec_max)
            max_max = np.max([ra_max, dec_max])

            if self.imagesize is None:
                self.imagesize = [max_max,max_max]

        # spatial_trim update begin

        if self.imagesize is None:
            max_size = np.max([len(self.xcoord),len(self.ycoord)])
            self.imagesize = [max_size,max_size]

        if np.array(self.imagesize).size == 1:
            self.imagesize=[self.imagesize,self.imagesize]

        wx,=np.where(np.abs(self.xcoord-self.obj_ra)*3600. * np.cos(np.deg2rad(self.obj_dec)) <= self.imagesize[0]) # edited by Eric
        wy,=np.where(np.abs(self.ycoord-self.obj_dec)*3600. <= self.imagesize[1])
        self.spatial_trim=[np.min(wx),np.max(wx)+1,np.min(wy),np.max(wy)+1]  # edited by Eric 

        # if self.imagesize is not None:
        #     if np.array(self.imagesize).size == 1:
        #         self.imagesize=[self.imagesize,self.imagesize]
        #     wx,=np.where(np.abs(self.xcoord-self.obj_ra)*3600. * np.cos(np.deg2rad(self.obj_dec)) <= self.imagesize[0]) # edited by Eric
        #     wy,=np.where(np.abs(self.ycoord-self.obj_dec)*3600. <= self.imagesize[1])
        #     self.spatial_trim=[np.min(wx),np.max(wx)+1,np.min(wy),np.max(wy)+1]  # edited by Eric 

        # else: # Edited by Eric
        #     if self.flat_cube is not None:
        #         mom0=(self.flat_cube > self.rmsfac*self.rms).sum(axis=2)
        #     else:
        #         mom0=(self.pbcorr_cube > self.rmsfac*self.rms).sum(axis=2)
        #     # mom0[mom0>0]=1
        #     # plt.figure();plt.imshow(mom0);plt.show() # Used by Eric to test
            
        #     cumulative_x = np.nancumsum(mom0.sum(axis=1),dtype=np.float)
        #     cumulative_x /= np.nanmax(cumulative_x)
        #     cumulative_y = np.nancumsum(mom0.sum(axis=0),dtype=np.float)
        #     cumulative_y /= np.nanmax(cumulative_y)
            
        #     wx_low,=np.where(cumulative_x < 0.02)
        #     if wx_low.size ==0: wx_low=np.array([0])
        #     wx_high,=np.where(cumulative_x > 0.98)
        #     if wx_high.size ==0: wx_high=np.array([cumulative_x.size])
        #     wy_low,=np.where(cumulative_y < 0.02)
        #     if wy_low.size ==0: wy_low=np.array([0])
        #     wy_high,=np.where(cumulative_y > 0.98)
        #     if wy_high.size ==0: wy_high=np.array([cumulative_y.size])
            
        #     beam_in_pix = int(np.ceil(self.bmaj/self.cellsize))

        #     # edited by Eric
        #     xpad = int(np.ceil(0.2 * (np.min(wx_high) - np.max(wx_low))))
        #     ypad = int(np.ceil(0.2 * (np.min(wy_high) - np.max(wy_low))))
        #     self.spatial_trim = [np.clip(np.max(wx_low) - 2*beam_in_pix - xpad,0,self.xcoord.size), np.clip(np.min(wx_high)+2*beam_in_pix+xpad,0,self.xcoord.size), np.clip(np.max(wy_low) - 2*beam_in_pix - ypad,0, self.ycoord.size), np.clip(np.min(wy_high) + 2*beam_in_pix + ypad,0,self.ycoord.size)]

        # spatial_trim end


        # check ra dec order?
        self.flat_cube_trim=self.flat_cube[self.spatial_trim[0]:self.spatial_trim[1],self.spatial_trim[2]:self.spatial_trim[3],self.chans2do[0]:self.chans2do[1]]
        self.pbcorr_cube_trim=self.pbcorr_cube[self.spatial_trim[0]:self.spatial_trim[1],self.spatial_trim[2]:self.spatial_trim[3],self.chans2do[0]:self.chans2do[1]]
        self.pb_trim = self.pb[self.spatial_trim[0]:self.spatial_trim[1], self.spatial_trim[2]:self.spatial_trim[3]] # added by Eric
        # self.mask_trim=self.mask[self.spatial_trim[0]:self.spatial_trim[1],self.spatial_trim[2]:self.spatial_trim[3],self.chans2do[0]:self.chans2do[1]] 
        self.spectralcube=self.spectralcube[self.chans2do[0]:self.chans2do[1],self.spatial_trim[2]:self.spatial_trim[3],self.spatial_trim[0]:self.spatial_trim[1]] 
        self.xcoord_trim=self.xcoord[self.spatial_trim[0]:self.spatial_trim[1]]
        self.ycoord_trim=self.ycoord[self.spatial_trim[2]:self.spatial_trim[3]]
        self.vcoord_trim=self.vcoord[self.chans2do[0]:self.chans2do[1]]  


    def smooth_mask(self,cube):  # updated by Eric
        """
        Apply a Gaussian blur, using sigma = 4 in the velocity direction (seems to work best), to the uncorrected cube.
        The mode 'nearest' seems to give the best results.
        :return: (ndarray) mask to apply to the un-clipped cube
        """
        sigma = self.spatial_smooth * self.bmaj / self.cellsize  # the version in use
        # sigma = self.spatial_smooth * np.sqrt(self.bmaj*self.bmin) / self.cellsize # another version considered

        cube[np.isnan(cube)] = 0.0 # the trimmed cube should be all within the ALMA FoV, thus not turning edge nan values to zero

        print(cube.shape,'cube shape before smoothing')
        # smooth_cube = ndimage.uniform_filter(cube, size=[sigma, sigma, self.spectral_smooth], mode='constant')  # mode='nearest' should be effectively the same, self.spectral_smooth was set to 4 channels
        smooth_cube = ndimage.gaussian_filter(cube,sigma=[sigma/2.35, sigma/2.35, self.spectral_smooth], mode='reflect')

        # noise_cube_big = smooth_cube[:,:, np.r_[0:self.chans2do[0], self.chans2do[1]:smooth_cube.shape[2]]] # very wrong
        noise_cube_big = self.flat_cube[:,:, np.r_[0:self.chans2do[0],self.chans2do[1]:self.flat_cube.shape[2]]]
        smooth_cube_noise=ndimage.gaussian_filter(noise_cube_big,sigma=[sigma/2.35, sigma/2.35,self.spectral_smooth],mode='constant')

        newrms = np.sqrt(np.nanmean(smooth_cube_noise**2))
        # newrms= self.rms_estimate(smooth_cube,0,1)  # commented by Eric; only use the first channel (the central square)
        print('RMS in smoothed cube (mJy/beam):',newrms*1000)

        self.cliplevel=newrms*self.rmsfac
        # self.cliplevel=self.rms*self.rmsfac

        mask=(smooth_cube > self.cliplevel)
        print('clip level (how many sigma)', self.rmsfac)

        if self.clean_mask is not None:
            hdu_clean = fits.open(self.clean_mask)
            clean_mask_data = np.squeeze(hdu_clean[0].data.T.astype(bool)) # transpose to match other cubes from read_in_a_cube()
            if self.flipped: clean_mask_data=np.flip(clean_mask_data,axis=2)
            clean_mask_data_trim = clean_mask_data[self.spatial_trim[0]:self.spatial_trim[1],self.spatial_trim[2]:self.spatial_trim[3],self.chans2do[0]:self.chans2do[1]]
            print('clean mask sum',np.sum(clean_mask_data_trim))
            mask = mask * clean_mask_data_trim

        mask_iter = mask.T # deliberately make them the same variable, convenient for updating
        
        if self.islandsize is None:
            self.islandsize = np.pi/(4*np.log(2)) * self.bmaj * self.bmin / self.cellsize**2 # beam size 
        if self.holesize is None:
            self.holesize = 2 * self.islandsize # (0.4 / self.cellsize)**2
        print('hole size threshold: %d spaxels, island size threshold: %d spaxels' %(self.holesize, self.islandsize))

        if self.holesize > 0 or self.islandsize>0:
            for i,v in enumerate(mask_iter):
                mask_iter[i] = morphology.remove_small_holes(mask_iter[i], self.holesize,connectivity=2)
                mask_iter[i] = morphology.remove_small_objects(v,self.islandsize,connectivity=2)

        mask_out = mask.astype(float).T

        index_dec = np.where(mask_out == True)[1]
        index_ra = np.where(mask_out == True)[2]
        ra_max = np.max(np.absolute( self.xc[index_ra]))
        dec_max = np.max(np.absolute( self.yc[index_dec]))
        print('Final mask, RA max (arcsec)', ra_max)
        print('Final mask, Dec max (arcsec)', dec_max)
        print('final mask sum',np.sum(mask))

        self.write_cube(mask_out)

        mask_out_full = np.pad(mask_out, pad_width=( (self.chans2do[0], self.flat_cube.shape[2]-self.chans2do[1])  
            , (self.spatial_trim[2], self.flat_cube.shape[1]-self.spatial_trim[3])
            , (self.spatial_trim[0], self.flat_cube.shape[0]-self.spatial_trim[1])))

        self.write_cube(mask_out_full, name='mask_full', version='no_trim')

        return mask      



    def make_all(self,pdf=False,fits=False):
        self.set_rc_params()

        # Edited by Eric
        # fig = plt.figure( figsize=(12,8))
        fig = plt.figure( figsize=(22,15))

        gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[2,1.1], hspace=0.0)
        # gs0.tight_layout(fig)  # edited by Eric
        gs00 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs0[0], wspace=0.0)

        ax1 = fig.add_subplot(gs00[0,0])
        ax2 = fig.add_subplot(gs00[0,1], sharey=ax1)
        ax3 = fig.add_subplot(gs00[0,2], sharey=ax1)
        textaxes = fig.add_subplot(gs00[0, 3])
        textaxes.axis('off')

        gs01 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1], wspace=0.3)

        ax4 = fig.add_subplot(gs01[0, 0])
        ax5 = fig.add_subplot(gs01[0, 1])

        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax3.get_yticklabels(), visible=False)
        # gs0.tight_layout(fig, rect=[0, 0, 1, 0.97]) # edited by Eric

        self.make_moments(axes=np.array([ax1,ax2,ax3]),fits=fits)
        self.make_pvd(axes=ax4,fits=fits)
        self.make_spec(axes=ax5,fits=fits)
        
        # ### plotting PA (rotation axis) on mom1
        # ypv=self.yc
        # xpv=self.yc*0.0
        # ang=360 - self.posang # edited by Eric to account for inversion of x-axis
        # c = np.cos(np.deg2rad(ang))
        # s = np.sin(np.deg2rad(ang))
        # x2 =  c*xpv - s*ypv
        # y2 =  s*xpv + c*ypv
        # ax2.plot(x2,y2,'k--')


        ###### make summary box
        
        rjustnum=35

        c = ICRS(self.obj_ra*u.degree, self.obj_dec*u.degree)

        thetext = (self.galname)+'\n \n'
        thetext += (("RA: "+c.ra.to_string(u.hour, sep=':',precision=2)))+'\n'
        thetext += ("Dec: "+c.dec.to_string(u.degree, sep=':', alwayssign=True,precision=2))+'\n \n'
        thetext += ("Vsys: "+str(int(self.vsys))+" km/s")+'\n'
        thetext += ("Dist: "+str(round(self.gal_distance,1))+" Mpc")+'\n'
        if self.gal_distance == self.vsys/70.:
            thetext+=("(Est. from Vsys)")



        at2 = AnchoredText(thetext,
                           loc='upper right', prop=dict(size=25,multialignment='right'), frameon=False, # textsize edited by Eric
                           bbox_transform=textaxes.transAxes
                           )
        textaxes.add_artist(at2)

        # plt.tight_layout()


        if pdf:
            plt.savefig(self.galname+"_allplots.pdf", bbox_inches = 'tight',dpi=50)
            plt.close()
        else:
            plt.show()


    def make_moments(self,mom=[0,1,2],axes=None,pdf=False,fits=False,circle=None,error=False,fov=None,mom0_angle=None,**kwargs):
        mom=np.array(mom)
        self.fits=fits
        
        if np.any(self.xc) == None:
            self.prepare_cubes()
        
        self.set_rc_params(mult=0.95) # edited by Eric
       
        nplots=mom.size
        if np.any(axes) == None:
            if self.make_square:
                fig,axes=plt.subplots(1,nplots,sharey=True,figsize=(7*nplots,12), gridspec_kw = {'wspace':0, 'hspace':0}) # size edited by Eric
            else:
                fig,axes=plt.subplots(1,nplots,sharey=True,figsize=(7*nplots*(self.imagesize[0]/self.imagesize[1]),7), gridspec_kw = {'wspace':0, 'hspace':0})
                
            outsideaxis=0
        else:
            outsideaxis=1
            
            
        if nplots == 1:
            axes=np.array([axes])

        store = []
        
        for i in range(0,nplots):
            if mom[i] == 0:
                mom_temp = self.mom0(axes[i],first=not i,last=(i==nplots-1),circle=circle,fov=fov,mom0_angle=mom0_angle)
                store.append(mom_temp)
            if mom[i] == 1:
                mom_temp = self.mom1(axes[i],first=not i,last=(i==nplots-1),error=error)
                store.append(mom_temp)
            if mom[i] == 2:
                mom_temp = self.mom2(axes[i],first=not i,last=(i==nplots-1),error=error)
                store.append(mom_temp)

            if mom[i] == 3:
                self.plot_mask(axes[i],first=not i,last=(i==nplots-1),**kwargs)
        
        # if self.gal_distance != None:
        #     self.scalebar(axes[i]) # plot onto last axes

        if pdf:
            if self.dpi is None:
                plt.savefig(self.galname+"_moment"+"".join(mom.astype(str))+".pdf", bbox_inches = 'tight') # 
                plt.close()
            else:
                plt.savefig(self.galname+"_moment"+"".join(mom.astype(str))+".png", bbox_inches = 'tight', dpi=self.dpi) # 
                plt.close()
        else:
            if not outsideaxis: plt.show()
        
        return store


    def plot_mask(self,ax1,first,last,**kwargs):

        try: channel = kwargs['channel']
        except: channel = 9

        try: title_text = kwargs['title_text']
        except: title_text = 'mask'

        im = ax1.imshow(self.mask_trim.T[channel], extent=[self.imagesize[0],-self.imagesize[0],-self.imagesize[1],self.imagesize[1]],interpolation='none')

        im.axes.tick_params(axis='both', which='both', color='white', top=True, right=True)
        
        axis_temp_x, axis_temp_y = self.axis_label(ax1, first, last)
        axis_temp_x.tick_params(axis='both', which='both', color='white', top=True, right=True)
        if last:
            axis_temp_y.tick_params(axis='both', which='both', color='white', top=True, right=True)

        # method 1 issue: can't find a way to hide the colorbar
        # cb=plt.gcf().colorbar(im, ax=[ax1],location='top',shrink=0.98,pad=0.09)  #         
        # cb.set_label(title_text, labelpad=10)

        # method 2 issue: either cannot match the positioning between colorbar() and divider or cannot shrink the bar in other panels with divider
        # cb_size = '%.5f%%' %( 100*(self.cb_pad + self.cb_frac)/(1-self.cb_pad-self.cb_frac) )
        # cb_pad = '%.5f%%' %(100* 0.0)
        # print(cb_size, cb_pad)
        # divider = make_axes_locatable(ax1)
        # cax = divider.append_axes("top", size=cb_size,pad=cb_pad)
        # plt.axis('off')

        # method 3 sucecess:
        axins1 = inset_axes(ax1, width="98%", height="5%", loc="upper center", bbox_to_anchor=(0., 0.22, 1.0, 1.0), bbox_transform=ax1.transAxes, borderpad=0,) # bbox_transform = (left, bottom, width, height)
        plt.axis('off')


        ax1.set_title(title_text,fontsize=50,pad=30) # normally 50, except for 35 for the third row

    def mom0(self,ax1,first=True,last=True,circle=None,fov=None,mom0_angle=None):

        mom0=(self.pbcorr_cube_trim*self.mask_trim).sum(axis=2)*self.dv
        print('With mask, total flux in Jy km/s',np.nansum(mom0/self.beam_area())) # added by Eric
        # mom0_raw=(self.pbcorr_cube_trim).sum(axis=2)*self.dv
        # print('Without mask, raw total sum flux in plotted region Jy km/s',np.nansum(mom0_raw/self.beam_area())) # added by Eric

        mom0[mom0<=0.] = np.nan
        mom0_err = self.rms / self.pb_trim * np.sqrt( self.mask_trim.sum(axis=2) ) * self.dv # added by Eric
        mom0_err[np.isnan(mom0)] = np.nan # added by Eric

        # testing code
        # mom0_err2 = np.zeros(self.pb_trim.shape)
        # for i,v1 in enumerate(self.pbcorr_cube_trim):
        #     for j,v2 in enumerate(v1):
        #         if mom0[i,j] > 0:
        #             spectrum_this = unumpy.uarray(self.pbcorr_cube_trim[i,j], self.rms/self.pb_trim[i,j])
        #             mom0_err2[i,j] = unumpy.std_devs( (spectrum_this * self.mask_trim[i,j]).sum() ) * self.dv
        # mom0_err2[mom0_err2<=0.] = np.nan


        oldcmp = cm.get_cmap("YlOrBr", 512)
        newcmp = ListedColormap(oldcmp(np.linspace(0.15, 1, 256)))

        # if self.mom0_tick is None: self.mom0_tick = np.linspace(0, (np.round((np.nanmax(mom0) / 10**np.floor(np.log10(np.nanmax(mom0)))))+1)*10**np.floor(np.log10(np.nanmax(mom0))), 5)
        # im1=ax1.contourf(self.xc, self.yc, mom0.T, levels=np.linspace(0, np.nanmax(mom0), 15), cmap=newcmp, zorder=1) # , norm=LogNorm(vmin=0.001,vmax=np.nanmax(mom0))  # tested by Eric np.max(self.mom0_tick) or self.maxmom0 , extend='max'

        im1 = ax1.imshow(mom0.T, vmin=0, vmax=self.maxmom0, cmap=newcmp, extent=[self.imagesize[0],-self.imagesize[0],-self.imagesize[1],self.imagesize[1]]) # 

        if mom0_angle is not None:
            if np.isscalar(mom0_angle):
                mom0_angle = np.array([mom0_angle])
            ### plotting PA (rotation axis) on mom1
            ypv=self.yc
            xpv=self.yc*0.0
            for v in mom0_angle:
                ang=180 - v # edited by Eric to account for inversion of x-axis
                c = np.cos(np.deg2rad(ang))
                s = np.sin(np.deg2rad(ang))
                x2 =  c*xpv - s*ypv
                y2 =  s*xpv + c*ypv
                if self.plot_axis:
                    ax1.plot(x2,y2,ls='--',color='magenta')

        # added by Eric
        if self.centre_marker:
            ax1.scatter(0,0,marker='x',color='cyan')

        
        # cb=plt.gcf().colorbar(im1, ax=[ax1],location='top',shrink=0.98,pad=self.cb_pad,fraction=self.cb_frac,ticks=self.mom0_tick)  # ,pad=0.09,fraction=0.15

        axins1 = inset_axes(ax1, width="98%", height="5%", loc="upper center", bbox_to_anchor=(0., 0.22, 1.0, 1.0), bbox_transform=ax1.transAxes, borderpad=0,) # bbox_transform = (left, bottom, width, height)
        cb = plt.gcf().colorbar(im1, cax=axins1, location='top',ticks=self.mom0_tick)

        # ax1_divider = make_axes_locatable(ax1)
        # cax2 = ax1_divider.append_axes("top", size=0.25, pad=1.0)
        # cax2_divider = make_axes_locatable(cax2)
        # cax2_divider.append_axes('left',size=0.1); cax2_divider.append_axes('right',size=0.1)
        # cb = plt.gcf().colorbar(im1, cax=cax2, location='top', ticks=self.mom0_tick)

        if self.bunit.lower() == "Jy/beam".lower(): cb.set_label("$I_{\\rm CO}$ (Jy beam$^{-1}$ km s$^{-1}$)", labelpad=10)
        if self.bunit.lower() == "K".lower(): cb.set_label("$I_{\\rm CO}$ (K km s$^{-1}$)", labelpad=10)
        
        self.add_beam(ax1)
        # self.aspect(ax1)
        self.axis_label(ax1, first, last)

        if circle is not None:
            if (type(circle)==int) or (type(circle) == float):
                circle = np.array([circle])
            for v in circle:
                radius = self.ang2pctrans_inv(v) #  / self.cellsize
                # circle1 = plt.Circle((0,0), radius, color='k',fill=False,linewidth=1, zorder=10)
                circle1 = Ellipse(xy=(0,0), width=radius*2, height=2*radius*np.cos(np.deg2rad(self.inc)), angle=(360-self.posang)-90, 
                    edgecolor='k', fc='None', lw=1,  zorder=10)
                ax1.add_patch(circle1)

        if fov is not None:
            if np.isscalar(fov):
                circle1 = Ellipse(xy=(0,0), width=fov[0], height=fov, angle=0, 
                    edgecolor='k', fc='None', lw=1, ls='-', zorder=10)
            else:
                all_colors = ['k','k','b','r']
                all_styles = ['-','--','-','-']
                for i in range(len(fov)):
                    circle1 = Ellipse(xy=(0,0), width=fov[i], height=fov[i], angle=0, 
                        edgecolor=all_colors[i], fc='None', lw=1, ls=all_styles[i], zorder=10)
                    ax1.add_artist(circle1)

        # plt.xlim(self.imagesize[0],-self.imagesize[0])
        # plt.ylim(-self.imagesize[0],self.imagesize[0])

        if self.fits:
            self.write_fits(mom0.T,0)
            self.write_fits(mom0_err.T,90) # added by Eric
            # self.write_fits(mom0_err2.T,80) # testing code

        return mom0.T # added by Eric
        

    def mom1(self,ax1,first=True,last=True,error=False):

        cube_temp = self.pbcorr_cube_trim.copy()
        # cube_temp[cube_temp<=0] = np.nan  # disagreed by Martin
        mom0=np.nansum(cube_temp*self.mask_trim, axis=2) # no need for dv here
        
        mom1 = np.nansum(cube_temp*self.mask_trim*self.vcoord_trim, axis=2) / mom0 # edited by Eric
        # mom1[mom1<=0] = np.nan; 
        mom1[~np.isfinite(mom1)] = np.nan
        
        if error:
            mom1_err2 = mom0.copy()*np.nan
            for i,v1 in enumerate(cube_temp):
                for j,v2 in enumerate(v1):
                    select = (~np.isnan(v2)) & (self.mask_trim[i,j] > 0)
                    selected = v2[select]
                    if len(selected) > 1: # (mom0[i,j] > 0.)
                        spectrum_this = unumpy.uarray(selected, self.rms/self.pb_trim[i,j])
                        mom1_err2[i,j] = unumpy.std_devs( np.sum(spectrum_this * self.vcoord_trim[select]) / np.sum(spectrum_this) )
            mom1_err2[mom1_err2==0] = np.nan; mom1_err2[~np.isfinite(mom1_err2)] = np.nan


        if self.max_velocity is None:
            self.max_velocity = np.nanmax(np.absolute(mom1.T-self.vsys))
        if self.range_spec is None:
            self.range_spec = self.max_velocity * 1.3 + 40

        # if self.vticks is None:
        #     self.vticks=np.linspace((-1)*np.ceil(np.max(np.abs(np.array(self.range_spec)-self.vsys))/10.)*10., np.ceil(np.max(np.abs(self.vcoord_trim-self.vsys))/10.)*10., 9)
        # vticks=np.linspace(self.range_spec[0]-self.vsys, self.range_spec[1]-self.vsys, 5)

        # levels=np.linspace(self.range_spec[0]-self.vsys, self.range_spec[1]-self.vsys, 500)
        # levels = np.linspace(-75, 75, 500)
        # im1=ax1.contourf(self.xc,self.yc, mom1.T-self.vsys, levels=levels, cmap=sauron ) # levels=self.vcoord_trim-self.vsys , vmin=self.vticks[0], vmax=self.vticks[-1]

        newcmp = ListedColormap(sauron(np.linspace(0.05, 0.95, 256)))
        
        im1 = ax1.imshow(mom1.T-self.vsys, vmin=-self.max_velocity, vmax=self.max_velocity, cmap=newcmp, extent=[self.imagesize[0],-self.imagesize[0],-self.imagesize[1],self.imagesize[1]]) # 

        print('Using vsys=%.2f km/s'%self.vsys)
        
        # edited by Eric
        # cb=self.colorbar(im1,ticks=vticks, pad=0.09, shrink=0.98)
        # cb=plt.gcf().colorbar(im1, ax=[ax1], location='top',shrink=0.98,pad=self.cb_pad, ticks=self.vticks)
        axins1 = inset_axes(ax1, width="98%", height="5%", loc="upper center", bbox_to_anchor=(0., 0.22, 1.0, 1.0), bbox_transform=ax1.transAxes, borderpad=0,) # bbox_transform = (left, bottom, width, height)
        cb = plt.gcf().colorbar(im1, cax=axins1, location='top', ticks=self.vticks)

        cb.set_label("$v_{\\rm "+ self.subscript +"} - v_{\\rm sys}$ (km s$^{-1}$)", labelpad=10) # edited by Eric

        
        self.add_beam(ax1)
        # self.aspect(ax1)
        self.axis_label(ax1, first, last)

        # added by Eric
        ### plotting PA (rotation axis) on mom1
        ypv=self.yc
        xpv=self.yc*0.0
        ang=180 - self.posang # edited by Eric to account for inversion of x-axis
        c = np.cos(np.deg2rad(ang))
        s = np.sin(np.deg2rad(ang))
        x2 =  c*xpv - s*ypv
        y2 =  s*xpv + c*ypv
        if self.plot_axis:
            ax1.plot(x2,y2,'k--')

        # added by Eric
        if self.centre_marker:
            ax1.scatter(0,0,marker='x',color='gray',zorder=100)
        
            
        if self.fits:
            self.write_fits(mom1.T,1)
            if error:
                self.write_fits(mom1_err2.T, 81)
                
        return mom1.T # added by Eric
        
    def mom2(self,ax1,first=True,last=True,error=False):

        cube_temp = self.pbcorr_cube_trim.copy()
        # cube_temp[cube_temp<=0] = np.nan # self.cliplevel  disagreed by Martin

        mom0=np.nansum(cube_temp*self.mask_trim, axis=2) # no need for dv here

        mom1=mom0.copy()*np.nan
        mom1=np.nansum(cube_temp*self.mask_trim*self.vcoord_trim, axis=2) / mom0 # edited by Eric
        mom1[~np.isfinite(mom1)] = np.nan

        # mom1[mom0 != 0.0] = (((self.pbcorr_cube_trim*self.mask_trim)*self.vcoord_trim).sum(axis=2))[mom0 != 0.0]/mom0[mom0 != 0.0]
        mom2=mom0.copy()*np.nan
        mom2_err2 = mom0.copy()*np.nan

        plot_count = 0
        for i in range(0,self.xcoord_trim.size):
            for j in range(0,self.ycoord_trim.size):
                select = (~np.isnan(cube_temp[i,j,:])) & (self.mask_trim[i,j,:] >0)
                selected = cube_temp[i,j,:][select]
                if (len(selected)>1) and (~np.isnan(mom1[i,j])):
                    number_this = len(select)
                    mom2[i,j]= np.sqrt(np.nansum( selected * (self.vcoord_trim[select] - mom1[i,j]) ** 2) / (np.nansum(selected) * ((number_this-1) / number_this ) ))

                    # The following block for testing
                    # print(np.nansum( selected * (self.vcoord_trim[select] - mom1[i,j]) ** 2) / (np.nansum(selected) * ((number_this-1) / number_this ) ), np.nansum(selected), np.nansum( selected * (self.vcoord_trim[select] - mom1[i,j]) ** 2))

                    # if np.nansum( selected * (self.vcoord_trim[select]-mom1[i,j]) ** 2) / (np.nansum(selected)*((number_this-1)/number_this)) < 0.: 
                    #     if (plot_count<10) and (np.nansum(selected)>0.03):
                    #         plt.figure()
                    #         plt.plot(selected,label='extracted spectrum within mask')
                    #         plt.legend()
                    #         plt.ylabel('flux')
                    #         plt.xlabel('channel')
                    #         ylim=np.max(np.absolute(selected))
                    #         plt.ylim(-ylim,ylim)
                    #         ax1 = plt.gca()
                    #         ax2 = ax1.twinx()
                    #         plt.plot(selected * (self.vcoord_trim[select]-mom1[i,j]) ** 2,label='spectrum multiplied with (v-vsys)^2',color='r')
                    #         plt.ylabel('flux * (v-vsys)^2')
                    #         plt.axhline(0,ls='--',color='gray')
                    #         ylim=np.max(np.absolute(selected* (self.vcoord_trim[select]-mom1[i,j]) ** 2))
                    #         plt.ylim(-ylim,ylim)
                    #         plt.legend(loc=2)
                    #         plt.show()
                    #         plot_count+=1
                    # if np.nansum(selected) < 0.:
                    # if np.nansum( selected * (self.vcoord_trim[select] - mom1[i,j]) ** 2) < 0.: 
                        # print('hhh updated')
                        # mom2[i,j] = 999.
                        # self.maxvdisp = 1005.
                    # Testing block ends.

                    if (error) and (~np.isnan(mom2[i,j])):
                        spectrum_this = unumpy.uarray(selected, self.rms/self.pb_trim[i,j])
                        vsys =  np.sum(spectrum_this * self.vcoord_trim[select]) / np.sum(spectrum_this)
                        mom2_err2[i,j] = unumpy.std_devs( unumpy.sqrt( np.sum( spectrum_this * (self.vcoord_trim[select] - vsys) ** 2) / (np.sum( spectrum_this) * ((number_this-1) /number_this ) ) ))

                # mom2[i,j] = len(selected) # for testing
                # self.maxvdisp = 5 # for testing; also comment out the next line 'mom2[mom2<=0.] ...'

        mom2[mom2<=0.] = np.nan; mom2[~np.isfinite(mom2)] = np.nan
        mom2_err2[mom2_err2<=0.] = np.nan; mom2_err2[~np.isfinite(mom2_err2)] = np.nan; mom2_err2[np.isnan(mom2)] = np.nan

        if self.maxvdisp == None:
            # self.maxvdisp = np.ceil(np.clip(np.nanstd(mom2)*4, 0, np.nanmax(mom2)) / 10.) * 10.
            self.maxvdisp=np.nanmax(mom2)
        
        # mom2levs=np.linspace(0,self.maxvdisp,15)
        # mom2levs2 = np.append(mom2levs, self.maxvdisp*100)
        # print(mom2levs2)
        
        # extend=None
        # if np.nanmax(mom2) > 1.2 * self.maxvdisp: extend='max'

        # im1=ax1.contourf(self.xc, self.yc, mom2.T, levels=mom2levs, cmap=sauron, vmax=self.maxvdisp*10) # , extend='max'

        newcmp = ListedColormap(sauron(np.linspace(0.05, 0.95, 256)))

        im1 = ax1.imshow(mom2.T, vmin=0, vmax=self.maxvdisp, cmap=newcmp, extent=[self.imagesize[0],-self.imagesize[0],-self.imagesize[1],self.imagesize[1]]) # 

        # added by Eric
        if self.centre_marker:
            ax1.scatter(0,0,marker='x',color='k')

        # vticks = np.linspace(0, (np.round((np.nanmax(mom2levs) / 10**np.floor(np.log10(np.nanmax(mom2levs))))) )* 10**np.floor(np.log10(np.nanmax(mom2levs))), 5)
        
        # edited by Eric
        # cb=self.colorbar(im1,ticks=vticks, pad=0.09, shrink=0.98)
        # cb=plt.gcf().colorbar(im1, ax=[ax1], location='top',shrink=0.98, pad=self.cb_pad)  # , ticks=vticks
        axins1 = inset_axes(ax1, width="98%", height="5%", loc="upper center", bbox_to_anchor=(0., 0.22, 1.0, 1.0), bbox_transform=ax1.transAxes, borderpad=0,) # bbox_transform = (left, bottom, width, height)
        cb = plt.gcf().colorbar(im1, cax=axins1, location='top')

        cb.set_label("$\sigma_{\\rm "+ self.subscript +"}$ (km s$^{-1}$)", labelpad=10)


        self.add_beam(ax1)
        # self.aspect(ax1)
        self.axis_label(ax1, first, last)

        if self.fits:
            self.write_fits(mom2.T,2)
            if error:
                self.write_fits(mom2_err2.T,82)

        return mom2.T # added by Eric

    def make_pvd(self,axes=None,fits=False,pdf=False,vel_file=None):
        
        self.fits = fits # added by Eric
        
        if np.any(self.xc) == None:
            self.prepare_cubes()
        
        if np.any(axes) == None:    
            self.set_rc_params(mult=0.95)   
            fig,axes=plt.subplots(1,figsize=(7,5)) # 8.5 height matching continuum image
            outsideaxis=0
        else:
            outsideaxis=1
                    
        # centpix_x=np.where(np.isclose(self.xc,0.0,atol=self.cellsize/1.9))[0]
        # centpix_y=np.where(np.isclose(self.yc,0.0,atol=self.cellsize/1.9))[0]

        cube_select = self.pbcorr_cube_trim*self.mask_trim

        # plt.close('all');plt.figure()
        # plt.imshow(cube_select[:,:,9], origin='upper'); plt.show()
        rotcube= rotateImage(cube_select,270-self.posang)  # ,[centpix_x[0],centpix_y[0]] # note that here the axes are [x,y,v], i.e. [x,y] are also reversed compared to common cases, thus 270-posang is correct
        # plt.close('all');plt.figure()
        # plt.imshow(rotcube[:,:,9], origin='upper'); plt.show()

        pvd=rotcube[:,np.array(rotcube.shape[1]//2-int(self.pvdthick/self.cellsize)).astype(int) : np.array(rotcube.shape[1]//2+int(self.pvdthick/self.cellsize)).astype(int),:].sum(axis=1)
        
        # changed by Eric
        loc1="upper left"
        # loc2="lower right"

        pvdaxis=(np.arange(0,pvd.shape[0])-pvd.shape[0]/2)*self.cellsize
        vaxis=self.vcoord_trim
        
        # pvd=pvd[np.abs(pvdaxis) < np.max([np.max(abs(self.xc)),np.max(abs(self.yc))]),:]
        # pvdaxis=pvdaxis[np.abs(pvdaxis) < np.max([np.max(abs(self.xc)),np.max(abs(self.yc))])]

        print('PVD shape', pvd.shape)
        
        oldcmp = cm.get_cmap("YlOrBr", 512)
        newcmp = ListedColormap(oldcmp(np.linspace(0.15, 1, 256)))
        
        axes.contourf(pvdaxis,vaxis,pvd.T,levels=np.linspace(self.cliplevel,np.nanmax(pvd),6),cmap=newcmp)
        axes.contour(pvdaxis,vaxis,pvd.T,levels=np.linspace(self.cliplevel,np.nanmax(pvd),6),colors='black')
        
        axes.set_xlabel('Offset (arcsec)')
        axes.set_ylabel('$v$ (km s$^{-1}$)')

        axes.tick_params(axis='y', which='both', right=False)
        secax = axes.secondary_yaxis('right', functions=(self.vsystrans, self.vsystrans_inv))
        secax.set_ylabel(r'$v - v_{\rm sys}$ (km s$^{-1}$)',rotation=270,labelpad=30)
        
        # added by Eric
        axes.tick_params(axis='x', which='both', top=False )

        if np.log10(self.ang2pctrans(np.max(np.concatenate([self.xc,self.yc])))) > 3: # edited by Eric
            secax2 = axes.secondary_xaxis('top', functions=(self.ang2kpctrans, self.ang2kpctrans_inv))
            secax2.set_xlabel(r'Offset (kpc)',labelpad=10)
        else:
            secax2 = axes.secondary_xaxis('top', functions=(self.ang2pctrans, self.ang2pctrans_inv))
            secax2.set_xlabel(r'Offset (pc)',labelpad=10)

        anchored_text = AnchoredText(r"$PA_{\rm kin} = $"+str(round(self.posang))+"$^{\circ}$", loc=loc1,frameon=False)
        axes.add_artist(anchored_text)

        # axes.set_ylim( axes.get_ylim()[0] * 0.995, axes.get_ylim()[1] * 1.005 )  # changed by Eric
        if self.range_spec is not None:
            axes.set_ylim( self.vsys - self.range_spec, self.vsys + self.range_spec)
        if self.pvd_radius is None:
            self.pvd_radius = self.imagesize[0]
        axes.set_xlim( -self.pvd_radius, self.pvd_radius) # added by Eric

        if vel_file is not None:
            radius_pc, vel_circ_rel = np.loadtxt(vel_file,comments=['#',';'],delimiter=',',unpack=True)
            radius_arcsec = self.ang2pctrans_inv(radius_pc)
            vel_circ_up = self.vsys + vel_circ_rel * np.sin(np.deg2rad(self.inc))
            vel_circ_low = self.vsys - vel_circ_rel * np.sin(np.deg2rad(self.inc)) 
            plt.plot(radius_arcsec,vel_circ_up,color='blue')
            plt.plot(-radius_arcsec,vel_circ_low,color='blue')

        # if self.gal_distance != None: # changed by Eric
        #     # self.scalebar(axes,loc=loc2)
        #     barlength_pc = 100.
        #     barlength_arc=  barlength_pc/(4.84*self.gal_distance)
        #     label="100 pc"
        #     asb = AnchoredSizeBar(axes.transData,  barlength_arc,   label,  loc="lower left",  pad=0.25, borderpad=0.5, sep=5, frameon=False)
        #     axes.add_artist(asb)

        if self.fits:
            self.write_pvd_fits(pvdaxis,vaxis,pvd.T)
            # print(pvd.T.shape)
        
        if pdf:
            plt.savefig(self.galname+"_pvd.pdf", bbox_inches = 'tight')
            plt.close()
        else:
            if not outsideaxis: plt.show()
    
    def make_spec(self,axes=None,fits=False,pdf=False,onlydata=False,onlymask=True,nsum=False,highlight=False,extra=False,label=None,update_v=False,high_color='orange',main_color='k',error_color='k',amp=1.0, mean_x=0.0, mean_y=0.0, fwhm=None, scaling=1.0, dialation=0):

        self.fits = fits # added by Eric

        if np.any(self.xc) == None:
            self.prepare_cubes()
        
        if axes == None:
            self.set_rc_params(mult=0.95)   
            fig,axes=plt.subplots(1,figsize=(8,6)) # (7.35,7.8) used for maser paper # (8,6) used for N1387 paper
            outsideaxis=0
        else:
            outsideaxis=1

        if extra == True:
            outsideaxis=1

        mask_spectrum = self.mask_trim.copy()
        if dialation>0:
            mask_spectrum = mask_spectrum.T
            for i in range(len(mask_spectrum)):
                mask_spectrum[i] = dia(mask_spectrum[i], iterations=dialation)
            mask_spectrum = mask_spectrum.T


        spec=self.pbcorr_cube_trim.sum(axis=0).sum(axis=0) * scaling
        spec_mask=(self.pbcorr_cube_trim*mask_spectrum).sum(axis=0).sum(axis=0) * scaling

        if fwhm is not None: # to synthesise single-dish observation
            size_x = self.pbcorr_cube_trim.shape[0]
            size_y = self.pbcorr_cube_trim.shape[1]
            size_v = self.pbcorr_cube_trim.shape[2]
            primary_beam = self.pri_beam(amp, mean_x, mean_y, fwhm, size_x, size_y)
            primary_beam_cube = np.tile(primary_beam, (size_v,1,1)).T # primary_beam (y,x) order, primary_beam_cube (x,y,v) order
            spec_mask=(self.pbcorr_cube_trim*mask_spectrum*primary_beam_cube * scaling).sum(axis=0).sum(axis=0)
            # print(primary_beam_cube.shape)
            # plt.figure()
            # plt.imshow(primary_beam)
            # plt.show()
            # plt.figure()
            # plt.imshow(primary_beam_cube[:,:,20])
            # plt.show()
            # plt.figure()
            # plt.imshow(mask_spectrum[:,:,20])
            # plt.show()

        # Noise calculation doesn't involve the scaling factor because that scaling only applies to emission lines, not to noise
        errors = np.zeros(len(spec))
        noise_cube = self.pbcorr_cube[ self.spatial_trim[0]:self.spatial_trim[1], self.spatial_trim[2]:self.spatial_trim[3], np.r_[0:self.chans2do[0],self.chans2do[1]:self.pbcorr_cube.shape[2]] ]
        print(noise_cube.shape[2],'channels used in error estimation')
        for i in range(len(errors)):
            # if np.nansum(mask_spectrum[:,:,i]) > 0.:
            if spec_mask[i] >0.:
                sum_all = np.zeros(noise_cube.shape[2])
                for j in range(len(sum_all)):
                    sum_all[j] = np.nansum(mask_spectrum[:,:,i] * noise_cube[:,:,j])
                errors[i] = np.nanstd(sum_all) # standard deviation
                # errors[i] = np.sqrt(np.sum(sum_all**2)/len(sum_all)) # root mean square
                # print('Given one line channel, mean of sums (Jy/beam):', np.nanmean(sum_all))
                # print('Given one line channel, std of sums (Jy/beam):', errors[i])
        # alternatively can use a sparse grid to make each spaxel independent and then calculate error propagation
        # print(errors)
        total_error = np.sqrt( np.sum( errors**2 ) )

        if self.bunit.lower() == "Jy/beam".lower():
            
            spec*=1/self.beam_area()
            spec_mask*=1/self.beam_area()
            errors*=1/self.beam_area()
            total_error*=1/self.beam_area()

            print('mean RMS per channel in synthesised spectrum (mJy):', np.mean(errors) * 1000)

            # added by Eric
            print(r'with mask, %.3f +- %.3f Jy km/s' %(np.sum(spec_mask)*self.dv, total_error*self.dv) )
            print('without mask, Jy km/s', np.sum(spec)*self.dv)
            # print(np.sum(mask_spectrum))
            # print(np.sum(mask_spectrum.T[10-self.chans2do[0]]))
            # plt.imshow(mask_spectrum.T[10-self.chans2do[0]]);plt.show()
            # plt.imshow(self.mask_trim.T[10-self.chans2do[0]]);plt.show()

            if spec_mask.max() < 1.0: # should be 1.0; temporarily increased to 30
                spec*=1e3
                spec_mask*=1e3
                errors*=1e3
                ylab="Flux density (mJy)"
            else:
                ylab="Flux density (Jy)"

        elif self.bunit.lower() == "k":
            ylab="Brightness Temp. (K)"
            print('with mask, K km/s',np.sum(spec_mask)*self.dv )
            print('without mask, K km/s', np.sum(spec)*self.dv)

        # added by Eric
        # left_v = self.vcoord_trim[ np.argmax( spec_mask[self.vcoord_trim < self.vsys])]
        # right_v = self.vcoord_trim[ np.argmax( spec_mask[self.vcoord_trim > self.vsys]) + np.sum(self.vcoord_trim <= self.vsys) ]
        # vsys = ( left_v + right_v ) / 2.

        vsys = np.nansum(spec_mask * self.vcoord_trim) / np.nansum(spec_mask)

        select = spec_mask>0.
        selected = spec_mask[select]
        spec_mask_uc = unumpy.uarray(selected, errors[select])
        vsys_err = unumpy.std_devs( np.sum(spec_mask_uc * self.vcoord_trim[select]) / np.sum(spec_mask_uc) )
        print('v_mean (mask applied): %.3f +- %.3f km/s' %(vsys,vsys_err) )

        # vsys_err_approx = np.sqrt(np.nansum(errors * self.vcoord_trim)**2) / np.nansum(spec_mask)
        # print('v error approx.', vsys_err_approx)

        if update_v == True:
            self.vsys = vsys

        if nsum:
            spec=np.append(running_mean(spec,nsum),spec[-1])
            spec_mask=np.append(running_mean(spec_mask,nsum),spec_mask[-1])
                
        # if self.range_spec is None:
        #     self.range_spec = [self.vcoord_trim[0]-20, self.vcoord_trim[-1]+20]
        # range_ind = (self.vcoord_trim >= self.range_spec[0])&(self.vcoord_trim <= self.range_spec[1])

        plt.xlim(self.vsys - self.range_spec, self.vsys + self.range_spec)

        if onlydata:
            axes.step(self.vcoord_trim,spec,c='k') # edited by Eric
            
            if highlight:
                plt.fill_between(self.vcoord_trim, spec, step="pre", alpha=0.3,color='grey')  # edited by Eric

        else:
            if not onlymask: # added by Eric
                axes.step(self.vcoord_trim,spec,c='lightgrey',linestyle='--') # edited by Eric

            if highlight:
                plt.fill_between(self.vcoord_trim,spec_mask, step="pre", alpha=0.3,color=high_color,label=label)
                axes.step(self.vcoord_trim,spec_mask,c=high_color)
                axes.plot([axes.get_xlim()[0],self.vcoord_trim[0]],[0,0],color=high_color)
                axes.plot([self.vcoord_trim[-1],axes.get_xlim()[1]],[0,0],color=high_color)
            else: # This is the default branch
                axes.step(self.vcoord_trim, spec_mask, c=main_color)
                axes.plot([axes.get_xlim()[0],self.vcoord_trim[0]],[0,0],color=main_color)
                axes.plot([self.vcoord_trim[-1],axes.get_xlim()[1]],[0,0],color=main_color)

            axes.fill_between(self.vcoord_trim, spec_mask-errors, spec_mask+errors, step='pre', color=error_color, alpha=0.5)


        axes.axhline(y=0,linestyle='dotted',color='k',alpha=0.5)        
        axes.set_xlabel('$v$ (km s$^{-1}$)') # _{\\rm LOS}
        axes.set_ylabel(ylab)
        
        axes.tick_params(axis='x', which='both', top=False )
        secax = axes.secondary_xaxis('top', functions=(self.vsystrans, self.vsystrans_inv))
        secax.set_xlabel('$v - v_{\\rm sys}$ (km s$^{-1}$)',labelpad=10) # edited by Eric
        
        if self.fits:
            self.write_spectrum(self.vcoord_trim,spec,self.vcoord_trim,spec_mask,ylab, errors)
        
        if pdf:
            plt.savefig(self.galname+"_spec.pdf", bbox_inches = 'tight')
            plt.close()
        else:
            if not outsideaxis: plt.show()
        
        return axes, self.vcoord_trim, spec_mask
        

    def make_continuum(self, cont_file, pdf=True, radius=None, vmax=None, sub=None, cont_angle=None, cont_beam=True): # added by Eric

        if np.any(self.xc) == None:
            self.prepare_cubes()

        self.set_rc_params(mult=0.95)

        continuum = fits.getdata(cont_file)
        hdr_cont = fits.getheader(cont_file,)
        unit = hdr_cont['BUNIT']
        if unit == 'Jy/beam':
            if np.nanmax(continuum) < 0.1:
                continuum *= 1000
                unit = r'(mJy beam$^{-1}$)'
            else:
                unit = r'(Jy beam$^{-1}$)'

        if vmax is None:
            vmax = np.nanmax(continuum)

        self.wcs=wcs.WCS(hdr_cont)
        print('\n check square pixel, check ra dec order: \n')
        print(self.wcs)

        xp = np.arange(hdr_cont['NAXIS1'])
        yp = np.ones(len(xp)) * int(hdr_cont['NAXIS2']/2)
        xcoord_cont = self.wcs.all_pix2world(xp,yp,0)[0]
        yp = np.arange(hdr_cont['NAXIS2'])
        xp = np.ones(len(yp)) * int(hdr_cont['NAXIS1']/2)
        ycoord_cont = self.wcs.all_pix2world(xp,yp,0)[1]

        select_ra = np.where(np.abs(xcoord_cont-self.obj_ra)*3600. * np.cos(np.deg2rad(self.obj_dec)) <= radius)[0]
        select_dec = np.where(np.abs(ycoord_cont-self.obj_dec)*3600. <= radius)[0]

        matplotlib.rc('axes',edgecolor='white')
        plt.figure(figsize=(7,12))

        im1 = plt.imshow(continuum[select_dec[0]:select_dec[-1]+1, select_ra[0]:select_ra[-1]+1], extent=[radius, -radius, -radius, radius], vmin=0, vmax=vmax) # 

        plt.gca().tick_params(axis='both', which='both', color='white', top=True, right=True)

        # if radius is not None:
        #     plt.xlim(-radius,radius);plt.ylim(-radius,radius)

        # plt.gca().set_aspect('equal') # equal is the default

        axis_temp_x, axis_temp_y = self.axis_label(plt.gca())
        axis_temp_x.tick_params(axis='both', which='both', color='white', top=True, right=True)
        axis_temp_y.tick_params(axis='both', which='both', color='white', top=True, right=True)

        if cont_beam:
            # self.add_beam(plt.gca(), edgecolor='white', facecolor='white')  # assuming continuum and line cube have the same beam
            bmaj_this=hdr_cont['BMAJ']*3600.
            bmin_this=hdr_cont['BMIN']*3600.
            bpa_this=hdr_cont['BPA']
            ae = AnchoredEllipse(plt.gca().transData, width=bmaj_this, height=bmin_this, angle=(360-bpa_this)-90, # corrected by Eric
                                 loc='lower left', pad=0.5, borderpad=0.4,
                                 frameon=False)                   
            ae.ellipse.set_edgecolor('white')
            ae.ellipse.set_facecolor('white')
            ae.ellipse.set_linewidth(1.5)
            plt.gca().add_artist(ae)

        if cont_angle is not None:
            if np.isscalar(cont_angle):
                cont_angle = np.array([cont_angle])
            ypv=np.linspace(-radius,radius,1000)
            xpv=ypv*0.0
            for v in cont_angle:
                ang=180 - v # edited by Eric to account for inversion of x-axis
                c = np.cos(np.deg2rad(ang))
                s = np.sin(np.deg2rad(ang))
                x2 =  c*xpv - s*ypv
                y2 =  s*xpv + c*ypv
                if self.plot_axis:
                    plt.gca().plot(x2,y2,ls='--',color='magenta')

        matplotlib.rc('axes',edgecolor='k')

        extend = None
        # if np.nanmax(continuum) > 1.2 * vmax:
        #     extend='max'

        cb=plt.gcf().colorbar(im1, ax=[plt.gca()], location='top',shrink=0.98,pad=self.cb_pad,extend=extend)
        if sub is None:
            cb.set_label(r'$I_{\rm continuum}$ '+unit, labelpad=10) # _{\rm'+freq+r'}
        else:
            cb.set_label(r'$I_{\rm '+sub+r'}$ '+unit, labelpad=10) # _{\rm'+freq+r'}


        if pdf==True:
            plt.savefig(self.galname+"_continuum.pdf", bbox_inches = 'tight')
            plt.close()
        else:
            plt.show()




    def scalebar(self,ax,loc='lower right'):  # edited by Eric
        # barlength_pc = np.ceil((np.abs(self.xc[-1]-self.xc[0])*4.84*self.gal_distance)/1000.)*100
        barlength_pc = np.ceil(1000.)
        barlength_arc=  barlength_pc/(4.84*self.gal_distance)
        
        if barlength_arc > 0.3*(self.xc[-1]-self.xc[0]): # go to 10 pc rounding
            # barlength_pc = np.ceil((np.abs(self.xc[-1]-self.xc[0])*4.84*self.gal_distance)/100.)*10
            barlength_pc = np.ceil(100.)
            barlength_arc=  barlength_pc/(4.84*self.gal_distance)

        if barlength_arc > 0.3*(self.xc[-1]-self.xc[0]): # go to 1 pc rounding
            # barlength_pc = np.ceil((np.abs(self.xc[-1]-self.xc[0])*4.84*self.gal_distance)/10.)*1
            barlength_pc = np.ceil(10.)
            barlength_arc=  barlength_pc/(4.84*self.gal_distance)
            
            
        
        if np.log10(barlength_pc) > 2.9: # slightly edited by Eric
            label=(barlength_pc/1e3).astype(str)+ " kpc"
        else:
            label=barlength_pc.astype(int).astype(str)+ " pc"
            
        asb = AnchoredSizeBar(ax.transData,  barlength_arc,   label,  loc=loc,  pad=0.25, borderpad=0.5, sep=5, frameon=False)
        ax.add_artist(asb)

    def add_beam(self,ax,edgecolor='black',facecolor='none'):
        ae = AnchoredEllipse(ax.transData, width=self.bmaj, height=self.bmin, angle=(360-self.bpa)-90, # corrected by Eric
                             loc='lower left', pad=0.5, borderpad=0.4,
                             frameon=False)                   
        ae.ellipse.set_edgecolor(edgecolor)
        ae.ellipse.set_facecolor(facecolor)
        ae.ellipse.set_linewidth(1.5)
        ax.add_artist(ae)

    def axis_label(self, ax1, first=True, last=True):

        ax1.tick_params(axis='y', which='both', right=False, left=False) 
        ax1.tick_params(axis='x', which='both', top=False )

        ax1.set_xlabel(r'RA - RA$\mathrm{_{centre}}$ (arcsec)') # 'RA offset (")' edited by Eric ,fontsize=self.lab_font
        if first:
            ax1.set_ylabel(r'Dec. - Dec.$\mathrm{_{centre}}$ (arcsec)') # 'Dec offset (")' ,fontsize=self.lab_font
            ax1.tick_params(axis='y', which='both', left=True)

        secax = None
        if np.log10(self.ang2pctrans(np.max(np.concatenate([self.xc,self.yc])))) > 3: # edited by Eric
            if last:
                secax = ax1.secondary_yaxis('right', functions=(self.ang2kpctrans, self.ang2kpctrans_inv))
                secax.set_ylabel(r'Dec. - Dec.$\mathrm{_{centre}}$ (kpc)',rotation=270,labelpad=30) # edited by Eric r'Dec offset (kpc)' ,fontsize=self.lab_font
            secax2 = ax1.secondary_xaxis('top', functions=(self.ang2kpctrans, self.ang2kpctrans_inv))
            secax2.set_xlabel(r'RA - RA$\mathrm{_{centre}}$ (kpc)',labelpad=10) # r'RA offset (pc)' ,fontsize=self.lab_font
        else:
            if last:
                secax = ax1.secondary_yaxis('right', functions=(self.ang2pctrans, self.ang2pctrans_inv))
                secax.set_ylabel(r'Dec. - Dec.$\mathrm{_{centre}}$ (pc)',rotation=270,labelpad=30) # edited by Eric r'Dec offset (kpc)' ,fontsize=self.lab_font
            secax2 = ax1.secondary_xaxis('top', functions=(self.ang2pctrans, self.ang2pctrans_inv))
            secax2.set_xlabel(r'RA - RA$\mathrm{_{centre}}$ (pc)',labelpad=10) # r'RA offset (pc)' ,fontsize=self.lab_font

        return secax2, secax

    def aspect(self,ax1):
        if self.make_square:
            largest_separation = np.max(np.absolute(np.concatenate([self.xc,self.yc])))
            # ax1.set_xlim(np.min([self.xc[0],self.yc[0]]),np.max([self.xc[-1],self.yc[-1]]))
            # ax1.set_ylim(np.min([self.xc[0],self.yc[0]]),np.max([self.xc[-1],self.yc[-1]]))
            ax1.set_xlim(largest_separation,-largest_separation)
            ax1.set_ylim(-largest_separation,largest_separation)
        ax1.set_aspect('equal')



    def write_cube(self, array, version='trim', name='mask'):
        filename = self.galname + '_'+name+'.fits'
        newhdu = fits.PrimaryHDU(array)

        newhdu.header['OBJECT']=self.spectralcube.header['OBJECT']
        try:
            newhdu.header['RADESYS']=self.spectralcube.header['RADESYS']
        except:
            pass
        try:
            newhdu.header['LONPOLE']=self.spectralcube.header['LONPOLE']
            newhdu.header['LATPOLE']=self.spectralcube.header['LATPOLE']
        except:
            pass
        try:
            newhdu.header['SPECSYS'] = self.spectralcube.header['SPECSYS']
        except:
            pass

        newhdu.header['VSYS']=(self.vsys,'km/s')
        newhdu.header['VELREF']=self.spectralcube.header['VELREF']
        newhdu.header['RESTFRQ']=self.spectralcube.header['RESTFRQ']

        pix1 = self.spectralcube.header['CRPIX1']
        # newhdu.header['CRPIX1']=pix1
        # newhdu.header['CRVAL1']=self.xcoord_trim[int(pix1-1)]  # self.spectralcube.header['CRVAL1']
        newhdu.header['CRPIX1']=1
        if version=='trim': newhdu.header['CRVAL1']=self.xcoord_trim[0]
        else: newhdu.header['CRVAL1']=self.xcoord[0]
        newhdu.header['CDELT1']=self.spectralcube.header['CDELT1']
        newhdu.header['CTYPE1']=self.spectralcube.header['CTYPE1']
        newhdu.header['CUNIT1']=self.spectralcube.header['CUNIT1']

        pix2 = self.spectralcube.header['CRPIX2']
        newhdu.header['CRPIX2']=1
        if version=='trim': newhdu.header['CRVAL2']=self.ycoord_trim[0]  # self.spectralcube.header['CRVAL1']
        else: newhdu.header['CRVAL2']=self.ycoord[0]
        newhdu.header['CDELT2']=self.spectralcube.header['CDELT2']
        newhdu.header['CTYPE2']=self.spectralcube.header['CTYPE2']
        newhdu.header['CUNIT2']=self.spectralcube.header['CUNIT2']

        pix3=self.spectralcube.header['CRPIX3']
        newhdu.header['CRPIX3']=1
        if version=='trim': newhdu.header['CRVAL3']=self.vcoord_trim[0]
        else: newhdu.header['CRVAL3']=self.vcoord[0]
        newhdu.header['CDELT3']=self.dv
        newhdu.header['CTYPE3']='VRAD'
        newhdu.header['CUNIT3']='km/s'

        newhdu.header['BMAJ']=(self.bmaj/3600.,'deg')
        newhdu.header['BMIN']=(self.bmin/3600.,'deg')
        newhdu.header['BPA']=(self.bpa,'deg')  # edited by Eric

        newhdu.header['spatial'] = (self.spatial_smooth, 'beam')
        newhdu.header['spectral'] = (self.spectral_smooth, 'channel')
        if np.isfinite(self.cliplevel): # added by Eric
            newhdu.header['RMSFAC']=(self.rmsfac)
            newhdu.header['MOMCLIP']=(self.cliplevel, self.bunit+' km/s')
        else:
            newhdu.header['RMSFAC']=('inf')
            newhdu.header['MOMCLIP']=('inf','inf')
        newhdu.header['hole'] = self.holesize
        newhdu.header['island'] = self.islandsize

        newhdu.header['comment'] = 'Created with pymakeplots'
        
        # print(newhdu.header)   
        newhdu.writeto(filename,overwrite=True)


    def write_fits(self,array,whichmoment):
        # added by Eric
        err = ''
        if int(str(whichmoment)[0]) > 5:
            err = '-err'+str(whichmoment)[0]
            whichmoment -= int(str(whichmoment)[0]) * 10

        # edited by Eric
        filename=self.galname+"_mom"+"".join(np.array([whichmoment]).astype(str))+err+".fits"
        
  
        newhdu = fits.PrimaryHDU(array)

        # changed by Eric
        newhdu.header['CRPIX1']=1
        newhdu.header['CRVAL1']=self.xcoord_trim[0]  # self.spectralcube.header['CRVAL1']
        newhdu.header['CDELT1']=self.spectralcube.header['CDELT1']
        newhdu.header['CTYPE1']=self.spectralcube.header['CTYPE1']
        newhdu.header['CUNIT1']=self.spectralcube.header['CUNIT1']
        newhdu.header['CRPIX2']=1
        newhdu.header['CRVAL2']=self.ycoord_trim[0] # self.spectralcube.header['CRVAL2']
        newhdu.header['CDELT2']=self.spectralcube.header['CDELT2']
        newhdu.header['CTYPE2']=self.spectralcube.header['CTYPE2']
        newhdu.header['CUNIT2']=self.spectralcube.header['CUNIT2']
        try:
            newhdu.header['PV2_1']=self.spectralcube.header['PV2_1']
            newhdu.header['PV2_2']=self.spectralcube.header['PV2_2']
        except:
            pass
            
        try:
            newhdu.header['RADESYS']=self.spectralcube.header['RADESYS']
        except:
            pass
            
        try:
            newhdu.header['SPECSYS'] = self.spectralcube.header['SPECSYS']
        except:
            pass    
        try:
            newhdu.header['LONPOLE']=self.spectralcube.header['LONPOLE']
            newhdu.header['LATPOLE']=self.spectralcube.header['LATPOLE']
        except:
            pass    
        newhdu.header['BMAJ']=self.bmaj/3600.
        newhdu.header['BMIN']=self.bmin/3600.
        newhdu.header['BPA']=self.bpa  # edited by Eric
        if np.isfinite(self.cliplevel): # added by Eric
            newhdu.header['MOMCLIP']=(self.cliplevel, self.bunit+' km/s')
        else:
            newhdu.header['MOMCLIP']=('inf','inf')
        newhdu.header['VSYS']=(self.vsys,'km/s')
        newhdu.header['comment'] = 'Moment map created with pymakeplots'

        if whichmoment == 0:
            newhdu.header['BUNIT']=self.bunit+' km/s'
        else:
            newhdu.header['BUNIT']='km/s'
                    
        newhdu.writeto(filename,overwrite=True)
    
    def write_pvd_fits(self,xx,vv,pvd):
        if self.fits == True:
            filename=self.galname+"_pvd.fits"
        else:
            filename=self.fits+"_pvd.fits"
            
        newhdu = fits.PrimaryHDU(pvd)
        newhdu.header['CRPIX1']=1
        newhdu.header['CRVAL1']=xx[0]
        newhdu.header['CDELT1']=xx[1]-xx[0]
        newhdu.header['CTYPE1']='OFFSET'
        newhdu.header['CUNIT1']='arcsec'
        newhdu.header['CRPIX2']=1
        newhdu.header['CRVAL2']=vv[0]
        newhdu.header['CDELT2']=vv[1]-vv[0]
        newhdu.header['CTYPE2']='VRAD'
        newhdu.header['CUNIT2']='km/s'
        newhdu.header['BMAJ']=self.bmaj/3600.
        newhdu.header['BMIN']=self.bmin/3600.
        newhdu.header['BPA']=self.bpa  # edited by Eric
        newhdu.header['PVDANGLE']=(self.posang,'deg')
        newhdu.header['PVDTHICK,half-width']=(self.pvdthick,'arcsec') # edited by Eric
        newhdu.header['MOMCLIP']=(self.cliplevel, self.bunit+' km/s')
        newhdu.header['VSYS']=(self.vsys,'km/s')
        newhdu.header['comment'] = 'Moment map created with pymakeplots'
        newhdu.header['BUNIT'] = self.bunit+' km/s'
        
        newhdu.writeto(filename,overwrite=True)               

    def write_spectrum(self,v1,spec1,vmask,specmask,descrip, errors):
       if self.fits == True:
           filename=self.galname
       else:
           filename=self.fits
       
       
       t = Table([v1,spec1],names=('Velocity (km/s)', descrip))   
       t.write(filename+"_spec.csv", format='csv',overwrite=True)
       
       t1 = Table([vmask,specmask, errors],names=('Velocity (km/s)', descrip, 'Error of flux'))   
       t1.write(filename+"_specmask.csv", format='csv',overwrite=True)




        # commented out by Eric


    # def colorbar(self,mappable,ticks=None,pad=0.05,**kwargs):
    #     ax = mappable.axes
    #     fig = ax.figure
    #     # divider = make_axes_locatable(ax)
    #     # cax = divider.append_axes("top", size="5%", pad=pad)
    #     # cb=fig.colorbar(mappable, cax=cax,ticks=ticks,orientation="horizontal",shrink=0.5)
    #     cb=fig.colorbar(mappable, ax=[ax], ticks=ticks, location='top',shrink=0.95,pad=pad)
    #     # cax.xaxis.set_ticks_position("top")
    #     # cax.xaxis.set_label_position('top')
    #     return cb        