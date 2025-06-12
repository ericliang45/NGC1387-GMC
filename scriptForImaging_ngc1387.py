import numpy as np
import os

from analysis_scripts import analysisUtils as au 
redshift = 0.00434
lambda_co = 299792458 / (230.538e9 / (1+redshift)) # in meters

# combined data
baseline = au.getBaselineLengths('NGC1387_combine.ms') # , field='NGC_3393'
lengths = np.array(np.array(baseline)[:,1], dtype=float)
print(np.median(lengths)) # 787.8180203595088
# Unprojected lengths:  min=8.902993, max=7552.388810, rms=1445.719265 m
l5 = lengths[ int(len(lengths) * 0.05) ]
lmin = np.min(lengths)
print(l5) # 5th percentile of uv distance (meters) L_5 = 71.8646
mrs1 = np.rad2deg(0.983 * lambda_co / l5) * 3600 # maximum recoverable scale in arcsec
print('MRS from 5th shortest baseline', mrs1) # 3.6848724
mrs2 = np.rad2deg(0.6 * lambda_co / lmin) * 3600
print('MRS from the shortest baseline', mrs2) # 18.1551271512


# ACA data
baseline = au.getBaselineLengths('uid___A002_Xc2bb44_X2142.ms.split.cal.split2.contsub') # , field='NGC_3393'
lengths = np.array(np.array(baseline)[:,1], dtype=float)
print(np.median(lengths)) # 21.10185995202084
# Unprojected lengths:  min=8.902993, max=44.725462, rms=23.834356 m
l5 = lengths[ int(len(lengths) * 0.05) ] 
lmin = np.min(lengths)
print(l5) # 5th percentile of uv distance (meters) L_5 = 9.085407069418887
mrs1 = np.rad2deg(0.983 * lambda_co / l5) * 3600 # maximum recoverable scale in arcsec
print('MRS from 5th shortest baseline', mrs1) # 29.146955772520833
mrs2 = np.rad2deg(0.6 * lambda_co / lmin) * 3600
print('MRS from the shortest baseline', mrs2) # 18.155127151206468



def getTargetSciSPW(thisMS):
    """ Function to get the science spectral windows (so ignoring #CH_AVG   """
    """ SPWs) for ALMA data.                                                """
    """ Returns science spectral window IDs as a comma separated string     """
    """ ready for use in e.g. split                                         """     
    
    """ To be run in CASA                                                   """
    """ Copy this function definition into your CASA script and run as      """
    """ mySPWstring =  getTargetSciSPW(myMS)                                """

    msmd.open(thisMS)
    targSPWs = msmd.spwsforintent('OBSERVE_TARGET#ON_SOURCE')
    chavgSPWs = msmd.almaspws(chavg=True,wvr=True)
    msmd.close()
    msmd.done()

    sciSPWstring = ','.join(map(str,np.setdiff1d(targSPWs,chavgSPWs)))

    return sciSPWstring


# parameters in step 2, average all time (1e9 sec), uv<100 (default meter), average all scans, for spw3 bin_chan=30 (18 km/s), for other spw, no chan bin.

myMS = 'uid___A002_Xc44eb5_X2e5d.ms.split.cal.split2'
contSPW_X2e5d = '0:0~70;90~127,1,2,3:0~1700;2200~3839' # X2e5d
# X2e5d, line at 231.3 GHz in spw 0, might be water v2=1 at (rest) 232.6867GHz

myMS = 'uid___A002_Xba26cb_X492.ms.split.cal.split2'
contSPW_X492 = '0,1:0~20;40~127,2,3:0~1700;2200~3839' # X492
# X492, line at 244.4 GHz in spw 1, rest frame 245.4 GHz

myMS = 'uid___A002_Xbbfdf3_X1224.ms.split.cal.split2'
contSPW_X1224 = '0,1:15~127,2,3:0~1700;2200~3839' # X1224

target = 'NGC1387'    #--- This is the target name as listed in the data.


myMS = 'uid___A002_Xc2bb44_X2142.ms.split.cal'
contSPW_ACA = '0:0~900;1150~2047,1:0~75;90~127,2,3'
target = 'NGC_1387'



nChanCS = 0                   #--- Number of channels to image in 13CS cube 
startCS = 0                    #--- Start Channel of CO cube . Format integer
myThreshCS = ''          #--- 13CS channel noise threshold. A string of format 'Y.ymJy

myChansRed = ''     #--- Redshifted emission channel range. A string of format 'X~Y'
myChansBlue = ''   #--- Redshifted emission channel range. A string of format 'X~Y'


#--- Logic switches ---#
""" Set the step you wish to run to = True and all others to = False """
step0 = False        #--- Run a listobs on the measurementSet myMS
step1 = False        #--- Find and split out science Spectral Windows (SPWs)
step2 = False        #--- Inpsect the data inPlotsMS
step3 = True        #--- Image the continuum.
step4 = False        #--- Image the CO and 13^CS line
step5 = False        #--- Moment maps



#-------------------#
#---- THE STEPS ----#

#-------------------------------------------------#
#--- STEP 0: Make listobs of the calibrated MS ---#
if step0:
    print "\n!!!STEP 0: Make listobs of the calibrated measurement set (MS).!!!\n"
    """ This step makes a listobs file of our calibrated MS. The listobs   """
    """ file contains useful information about the observation setup e.g.  """ 
    """ frequency settings, time on source, number of antennas and more.   """
    """ Making a listob at this stages means we can see what data in this  """
    """ MS we wish to retain.                                              """
    
    if os.path.exists(myMS):
        print '>> Generating listobs file of '+myMS
        listobs(myMS,listfile=myMS+'.listobs')

        """>> You can now inspect the listobs file in a standard text editor. """
    else:
        print ">> You do not seem to have "+myMS+" in the current directory, "
        print ">> skipping the rest of this step. """

#---------------------------------------------------------#
#--- STEP 1: Split out science target and science SPWs ---#
if step1:
    print "\n!!! STEP 1: Split out science target and science Spectral Windows (SPWs).!!!\n"
   
    """ myMS contains lots of data we won't actually need during imaging the"""
    """ science target. For example scans on the observed calibrators and   """
    """ spectral windows (SPWs) which are used only for calibration.        """
    """ To give us a smaller MS to work with we will first split out from   """
    """ myMS, the Target source and the science SPWs.                       """

    scispws = getTargetSciSPW(myMS)
    print '>> scispws: '+scispws

    split(vis=myMS,
        intent = 'OBSERVE_TARGET#ON_SOURCE',
        datacolumn='data', # corrected
        field=target,
        spw = scispws,
        outputvis=myMS+'.split2',)

    print '>> New measurement set: '+target+'_TM2.split.cal'
    print '>> Generating listobs file of new MS'

    """ Now we generate a listobs file to see the contents of our new       """
    """ measurement set                                                     """

    listobs(target+'_TM2.split.cal',listfile=target+'_TM2.split.cal.listobs')


#------------------------------------------#
#--- STEP 2: Inspect the data in plotMS ---#
if step2:
    print "\n!!! STEP 2: Inspect the data in plotMS !!!\n"
    print ">> This is a hands-on step, please see the instructions on the tutorial webpage."
   
    print ">> Once you have inspected the data please fill in the contSPW parameter in"
    print ">> the 'Initial parameters' section of this script."


#-------------------------------#
#--- STEP 3: Image continuum ---#
if step3:
    print "\n!!!STEP 3: Image continuum!!!\n"
    print ">> Here we will image the spectral line free continuum emission of "
    print ">> the target source Z-CMa.\n"
    print ">> Please read STEP3 of the tutorial, calculate and input the required"
    print ">> parameters into the 'Initial parameters' section of this script.\n"
    print ">> Before proceeding, please confirm you have set contSPW, myCell,"
    print ">> myImsize and myThreshold? [y/n]"

    # for X2e5d, Nate = 48, Ton = 725.5s, Band width = 7826170000 Hz, Tsys = 80K; X1224 100K, X492 90K
    # theoretical rms = 2 * 1.38064852e-23 * 80 / (np.pi * 6**2 * 0.7) / np.sqrt(48 * 47 * 725.5 * 7826170000) * 1e26 * 1e3 =0.0247 mJy
    # measured rms
    # myThreshold = '0.075mJy' #--- 3x your calculated theoretical threshold. A string of format 'Y.ymJy'

    # 1, threshold='0.09mJy'
    # 2, threshold='0.03mJy', update mask path
    tclean(
      vis=['tracks/uid___A002_Xba26cb_X492.ms.split.cal.split2',
    'tracks/uid___A002_Xbbfdf3_X1224.ms.split.cal.split2',
    'tracks/uid___A002_Xc44eb5_X2e5d.ms.split.cal.split2',
    'tracks/uid___A002_Xc2bb44_X2142.ms.split.cal.split2',],
      imagename = 'NGC1387_continuum_2',
      spw=[contSPW_X492, contSPW_X1224, contSPW_X2e5d, contSPW_ACA],
      threshold='0.03mJy', # measured rms = 0.02
      niter=1000000,
      cell=['0.04arcsec'],
      imsize=[750,750], #should cover about the FoV 
      outframe='bary',
      restfreq = '230.538GHz',
      deconvolver='multiscale',
      scales=[0,4,12],# inpixels equivalent to point source and 1, 3 times the beam.
      weighting = 'briggs',
      robust=1.0,
      pbcor=True, 
      specmode='mfs',
      restoringbeam='common',
      nterms=1,        
      chanchunks=-1, 
      gridder='mosaic',   
      interactive=False,
      mask = 'jacob_continuum/NGC1387_lower_continuum.clean.mask',
        )

rootfile = '/Users/liangf/work/ngc1387/calibrated_all/continuum_trial/NGC1387_continuum_2.image.pbcor'
exportfits(
    imagename=rootfile,
    fitsimage=rootfile+'.fits',
    dropstokes=True,
    dropdeg=True,
    )


#----------------------------------------------------#
#--- STEP 4: Subtract the continuum & Image lines ---#
if step4:
    print "\n!!!STEP 4: Subtract the continuum & Image lines!!!\n"        

    """ First we need to generate a continuum subtracted measurement set.   """
    """ For this we use the same SPW & channel ranges from the continuum    """
    """ imaging to give CASA a range of channels to fit the continuum to    """
    """ within the task uvcontsub"""

    myMS = 'uid___A002_Xc2bb44_X2142.ms.split.cal.split2'

    print ">> Generating continuum subtracted measurement set."
    uvcontsub(
        # vis = target+'_TM2.split.cal',
        vis = myMS,
        fitspw = contSPW,
        combine = 'spw',
        # field = 'NGC1387',
        spw='0',)
    print ">> "+myMS+'.contsub now created.'
    # default values: fitorder = 0, solit = 'int'. see document of uvcontsub_old https://casadocs.readthedocs.io/en/latest/api/tt/casatasks.manipulation.uvcontsub_old.html

    """ This creates a new MS called target+'_TM2.split.cal.contsub which we"""
    """ should use for line imaging                                         """

    """ Now we image the CO line, we will use CASA's automasking capability."""


    """ CLEANING OF THE CO CUBE                                             """
    #---uncommeny below lines if you want to delete existing image files ---#
    #--- with same file name before imaging  ---#
    #os.system('rm -rf  '+target+'_CO_cube.image')
    #os.system('rm -rf  '+target+'_CO_cube.image.pbcor')
    #os.system('rm -rf  '+target+'_CO_cube.psf')
    #os.system('rm -rf  '+target+'_CO_cube.mask')
    #os.system('rm -rf  '+target+'_CO_cube.model')
    #os.system('rm -rf  '+target+'_CO_cube.residual')
    #os.system('rm -rf  '+target+'_CO_cube.sumwt')
    #os.system('rm -rf  '+target+'_CO_cube.pb')

    myMS1 = 'uid___A002_Xc44eb5_X2e5d.ms.split.cal.split2'
    myMS2 = 'uid___A002_Xba26cb_X492.ms.split.cal.split2'
    myMS3 = 'uid___A002_Xbbfdf3_X1224.ms.split.cal.split2'
    target = 'NGC1387'    #--- This is the target name as listed in the data.

    '''
    import analysisUtils as au
    blines = au.getBaselineLengths(myMS1)
    '''
    # FOV = 299792458 / 230.5380e9 / 12. / np.pi * 180 * 3600 = 22.35 arcsec (24.44" on ALMA archive, FWHM=27.4" in manual)
    # Resolution = 299792458 / 230.5380e9 / 7552.3888096005703 / np.pi * 180 * 3600 = 0.0355 arcsec for X2e5d （0.047" on ALMA archive）, m/s/Hz/m
    # Sell size = resolution / 5 = 0.0355 / 5 = 0.007 arcsec
    # Sell size = resolution / 3 = 0.0355 / 3 = 0.012 arcsec
    # Imsize = FOV / Sell size * 2 = 22.35 / 0.007 * 2 = 6386
    # Read from available cubes, FoV of the clouds are within 20" (DEC) * 25.5" (RA), 26 / 0.012 = 2167, 26/0.007=3714~3840=256*3*5
    # The requirement is even number and factorizable by 2,3,5,7 only. 6386 ~  4 * 9 * 25 *  7 = 6300
    myCell = '0.007arcsec'  #--- cell or spaxel size of your images. A string of format 'X.xxxarcsec' 
    myImsize = 3840  #--- number of x by y pixels in your image. Format integer 

    startCO = '1170km/s' #                   #--- Start Channel of CO cube. Format integer
    # 229.5 - 229.65 GHz
    # Radio velocity
    # v = c*(1-f/f0)=299792.458*(1- (230580.444e3 + 900*(-976.562)) / 230.5380e6) = 1087 km/s
    # v = c*(1-f/f0)=299792.458*(1- (230580.444e3 + 1150*(-976.562)) / 230.5380e6) = 1405 km/s
    # Read from high-res cube, line exist from 1030km/s to 1210km/s

    nChanCO = 100 #--- Number of channels to image in CO cube. Format integer. 
    # 0.6 km/s * (2300 - 1600) / 2 km/s = 210, (1520.554 - 1072.228) / 2 = 224
    # Number of channels in the output image

    # 0.075 * np.sqrt(7826170000) / np.sqrt(488.281)
    myThreshCO = '300mJy'          #--- CO channel noise threshold. A string of format 'Y.ymJy'

    print "  >> Cleaning CO emission in Z-CMa."
    tclean(
        vis = [myMS1+'.contsub',myMS2+'.contsub',myMS3+'.contsub'],
        imagename = target+'_CO_cube',
        # field = target,
        stokes = 'I',
        # spw = '3', 
        outframe = 'LSRK',
        restfreq = '230.538GHz',
        specmode = 'cube',
        nchan = nChanCO,
        start = startCO, 

        width = '2km/s',
        cell=[myCell],
        imsize=[myImsize,myImsize], #covers about the FoV 
        gridder = 'standard',
        pbcor = True,
        restoringbeam = '',

        weighting = 'briggs',
        robust = 0.5,

        deconvolver = 'multiscale',
        scales=[0,5,15],# inpixels equivalent to point source and 1, 3 times the beam.
        niter = 0, # 10000000, 
        threshold = myThreshCO, 
        interactive = False,

        usemask = 'auto-multithresh',
        sidelobethreshold = 2.0,
        noisethreshold = 4.25,  
        minbeamfrac =  0.3,           
        lownoisethreshold = 1.5,
        negativethreshold = 0.0,

        chanchunks=-1,
        )

    ''' For ACA data alone '''

    # FOV = 299792458. / 230.5380e9 / 7. / np.pi * 180 * 3600 = 38.32 arcsec (41.847" on ALMA archive, FWHM=46.65 in manual)
    # Resolution = 299792458 / 230.5380e9 / 44.725462 / np.pi * 180 * 3600 = 5.9 arcsec（4.9" on ALMA archive）, m/s/Hz/m
    # Sell size = resolution / 5 = 5.9 / 6. = 1 arcsec
    # Imsize = FOV / Sell size * 2 = 38.32 / 1. * 2 = 76.6 ~ 81 
    # The requirement is to use multiples of 2, 3, 5, and 7
    myCell = '1.0arcsec'  #--- cell or spaxel size of your images. A string of format 'X.xxxarcsec' 
    myImsize = 84  #--- number of x by y pixels in your image. Format integer 

    startCO = '1170km/s' #                   #--- Start Channel of CO cube. Format integer
    # Optical velocity, useless
    # (230.5380e6 / (230497.653e3 + 1600 * (-488.281)) -1)* 299792.458 = 1072.228 km/s
    # (230.5380e6 / (230497.653e3 + 2300 * (-488.281)) -1)* 299792.458 = 1520.554 km/s
    # Radio velocity
    # v = c*(1-f/f0)=299792.458*(1- (230497.653e3 + 1600*(-488.281)) / 230.5380e6) = 1068 km/s
    # v = c*(1-f/f0)=299792.458*(1- (230497.653e3 + 2300*(-488.281)) / 230.5380e6) = 1512 km/s
    # Read from high-res cube, line exist from 1030km/s to 1210km/s

    nChanCO = 100 #--- Number of channels to image in CO cube. Format integer. 
    # Number of channels in the output image

    # ACA, dirty cube
    # ACA_clean[0-2], roughly the same cleaning. threshold=66mJy
    # ACA_clean3, threshold=21mJy, used old mask
    # ACA_clean4, threshold=14mJy, used auto-generated mask (masks are roughly the same)
    # ACA_clean5, large cube size, 2.5*. thre=14mJy. This doesn't work. Everything outside the primary beam is blank. The above cubes already touched the edge.

    # for ACA, Nate = 9, Ton = 1287.98s, Band width = 976562 Hz, Tsys = 65 K
    # theoretical rms = 2 * 1.38064852e-23 * 65 / (np.pi * 3.5**2 * 0.7) / np.sqrt(9 * 8 * 1287.98 * 976562) * 1e26 * 1e3 = 22.1 mJy， *3 recommended in the workshop
    # 14.83 mJy/beam from measurement of line free channels of dirty cube, *1.5 suggested by Martin
    myThreshCO = '14mJy'          # A string of format 'Y.ymJy'
    myMS = 'uid___A002_Xc2bb44_X2142.ms.split.cal.split2'
    target = 'NGC_1387'    #--- This is the target name as listed in the data.

    print "  >> Cleaning CO emission in Z-CMa."
    tclean(
        vis = myMS+'.contsub',
        imagename = target+'_cube_ACA_clean5',
        # field = target,
        stokes = 'I',
        # spw = '3', 
        outframe = 'LSRK',
        restfreq = '230.538GHz',
        specmode = 'cube',
        nchan = nChanCO,
        start = startCO, 

        width = '2km/s',
        cell=[myCell],
        imsize=[myImsize,myImsize], #covers about the FoV 
        gridder = 'standard',
        pbcor = True,
        restoringbeam = '',

        weighting = 'briggs',
        robust = 0.5,

        deconvolver = 'multiscale',
        scales=[0,5,15],# inpixels equivalent to point source and 1, 3 times the beam.
        niter = 10000000, # 10000000, 
        threshold = myThreshCO, 
        interactive = True,

        # mask = 'NGC_1387_CO_cube_ACA_clean2.mask', usemask='user',
        usemask = 'auto-multithresh',
        sidelobethreshold = 1.35,
        noisethreshold = 4.5,  # 5.0 recommended
        minbeamfrac =  0.1,           
        lownoisethreshold = 1.8, # 2.0 recommended
        negativethreshold = 0.0,
        fastnoise=True,

        # chanchunks=-1,
        )


#-------------------------------#
#--- STEP 5: Moment Analysis ---#
if step5:
    print '\n!!! STEP 5: Moment Analysis !!!\n'

    """ Finally we generate 0th, 1st and 2nd order image moments of the CO  """
    """ emission from Z-CMa. This gives us integrated line intensity,       """
    """ the velocity fields and the velocity distributions of the CO.       """

    print "  >> Now Generating redshifted CO moments 0,1 & 2 for Z-CMa."    
    immoments('NGC_1387_CO_cube_ACA.image.pbcor',
        moments =[0],
        axis='spectral',
        chans = '9~86', # 9~86 for ACA NGC_1387_CO_cube_ACA_clean0.image.pbcor 11~88 for 12-m cube  NGC1387_CO_cube.image.pbcor
        # excludepix = [-100000.0,0.0],
        outfile = 'NGC_1387_CO_cube_ACA.image.pbcor.moment'
        )

    # print "  >> Now Generating blueshifted CO moments 0,1 & 2 for Z-CMa."
    # immoments(imagename  = 'Z_CMa_CO_cube.image',
    #     moments =[0,1,2],
    #     axis='spectral',
    #     chans =myChansBlue,
    #     excludepix = [-100000.0,0.2],
    #     outfile = 'Z_CMa_CO_cube.image_blue'
    #     )
#--- END THE STEPS ---#
#---------------------#


exportfits(
    imagename='NGC_1387_cube_ACA_clean4.image.pbcor',
    fitsimage='NGC_1387_cube_ACA_clean4.image.pbcor.fits',
    dropstokes=True,
    dropdeg=True,
    )


''' check self-calibration '''

# 'uid___A002_Xba26cb_X492.ms.split.cal.split2', contSPW = '0,1:0~20;40~127,2,3:0~1700;2200~3839'
#     'uid___A002_Xbbfdf3_X1224.ms.split.cal.split2', contSPW = '0,1:15~127,2,3:0~1700;2200~3839'
#     'uid___A002_Xc44eb5_X2e5d.ms.split.cal.split2', contSPW = '0:0~70;90~127,1,2,3:0~1700;2200~3839'
#     'uid___A002_Xc2bb44_X2142.ms.split.cal.split2', contSPW = '0:0~900;1150~2047,1:0~75;90~127,2,3'

myCell = '0.1arcsec'
myImsize = 300

contSPW = '0,1:0~20;40~127,2,3:0~1700;2200~3839'

tclean(
    vis = 'uid___A002_Xba26cb_X492.ms.split.cal.split2',
    imagename = 'X492.continuum.clean2',
    spw=contSPW,
    cell=[myCell],
    imsize=[myImsize,myImsize], #should cover about the FoV 
    outframe='LSRK',
    deconvolver='multiscale',
    scales=[0,5,15],# inpixels equivalent to point source and 1, 3 times the beam.
    weighting = 'briggs',
    robust=1.0,
    pbcor=False, 
    specmode='mfs',
    restoringbeam='common',
    nterms=1,        
    chanchunks=-1, 
    gridder='standard',   
    interactive=True,

    niter=10000, # 10000
    threshold='60mJy',
    usemask = 'auto-multithresh',
    pbmask = 0.2,
    sidelobethreshold = 1.25,
    noisethreshold = 5.0,
    minbeamfrac =  0.1,  
    lownoisethreshold = 2.0,
    negativethreshold = 0.0,
    fastnoise=False,

    )

''' check self-calibration (end) '''


''' Combine 7-m and 12-m data '''

plotms(vis='',yaxis='wt',xaxis='uvdist',plotfile='XXX_WT.png',showgui=True)

concat(vis=['uid___A002_Xba26cb_X492.ms.split.cal.split2.contsub',
    'uid___A002_Xbbfdf3_X1224.ms.split.cal.split2.contsub',
    'uid___A002_Xc44eb5_X2e5d.ms.split.cal.split2.contsub',
    'uid___A002_Xc2bb44_X2142.ms.split.cal.split2.contsub',],
       concatvis='NGC1387_combine.ms')

listobs('NGC1387_combine.ms',listfile='NGC1387_combine.ms.listobs')

plotms(vis='NGC1387_combine.ms',yaxis='wt',xaxis='uvdist',coloraxis='spw',plotfile='combine_CO_WT.png',showgui=True)

plotms(vis='NGC1387_combine.ms',yaxis='amp',xaxis='uvdist',spw='', avgscan=True,
       avgchannel='5000', coloraxis='spw',plotfile='combine_CO_amp.png',showgui=True)


import analysisUtils as au
blines = au.getBaselineLengths('NGC1387_combine.ms')

# Synthesized beam. 0.0355" for longest baseline. 
# 15 pc resolution <-> 0.1738", sampled by 5 pixels, each pixel = 0.0348" ~ 0.035"; sampled by 3 pixels, 0.058"
# "NGC1387_CO_cube.image"(no ACA), robust=-0.5, 0.0923148, 0.0653339", 1.279 mJy
# "NGC1387_combine_dirty", "clean0"，robust=1.5, 0.223116, 0.174808, -63.2071 deg, 1.156 mJy
# "dirty1", small size, robust=0.5, 0.143167, 0.12516, 1.228 mJy # slightly dependent on cube size
# "dirty3", normal size, robust=0.5, 0.131743, 0.109433, 86.8341 deg, 1.171 mJy
# "dirty2", robust=0.5, uvtaper=0.1", 0.213733, 0.192974, 89.8488 deg, 1.223 mJy
# "dirty4", robust=0.5, uvtaper=0.18", 0.326244, 0.292959, 86.138 deg, 1.381 mJy
# "dirty5", robust=0.5, uvtaper=0.14", 0.282451" X 0.23596", 87.8165 deg, 1.300 mJy
# "dirty6", robust=-0.5, 0.103724" X 0.0553723", 50.2258 deg, 1.7380 mJy
# "dirty7", robust=0, 0.0963973" X 0.064646", 54.1452 deg, 1.40423 mJy
# "dirty8", robust=0 (not plotted), taper=0.22, 0.328636" X 0.299003", 67.9565 deg, 1.61378083184 mJy
# "dirty9", robust=0.5, uvtaper=0.22", 0.369997" X 0.333697", 84.0167 deg, 1.45590876254 mJy
# "dirty10", robust=1.5,  0.192152" X 0.148708", -74.9322 deg, 1.10353626367 mJy
# "dirty11", robust=-2 (not plotted), 0.099079" X 0.0488067", 56.5701 deg, 2.00337813961
# "dirty12", robust=2, 0.193312" X 0.149507", -74.7668 deg, 1.10353157787
# "dirty13", robust=1, 0.167105" X 0.135577", 89.5616 deg, 1.10739086562
# "dirty14", robust=1, uvrange='>200'(m), 0.137108" X 0.0940404", 79.4782 deg, 1.22756493811
# "dirty15", robust=1, uvrange='>100', 0.147594" X 0.127897", -84.6488 deg, 1.14780750372
# "dirty16", robust=1, uvrange='>300', 0.105005" X 0.0891811", -89.8258 deg, 1.3088422984
# "dirty17", robust=1, uvrange='<5000', 0.167105" X 0.135577", 89.5616 deg, 1.10739086562
# "dirty18", robust=1, uvrange='<4000', 0.167461" X 0.135786", 89.588 deg, 1.10834843061, 
# "dirty19", robust=1, uvrange='<3000', 0.18988" X 0.150236", -77.5519 deg, 1.14270023029
# "dirty20", robust=1, uvrange='<2000', 0.21855" X 0.171596", 86.7588 deg, 1.19528804916
# "dirty21", robust=1, uvrange='<1000', 0.343776" X 0.290804", -86.3688 deg, 1.4079168572
# "dirty22", robust=1, uvrange='<1250', 0.298649" X 0.244991", -85.4948 deg, 1.33868466959
# "dirty23", robust=1, uvrange='<1500', 0.261209" X 0.208502", -84.2532 deg, 1.26673033577
# "dirty24", robust=1, uvrange='<1750', 0.24185" X 0.19343", -86.1266 deg, 1.2238299519


myCell_dirty = '0.04arcsec'
myImsize_dirty = 750
# 30" / 0.058" = 512

robust=1

nChanCO = 100
startCO = '1170km/s'

# dirty cube
tclean( vis = 'NGC1387_combine.ms',
        imagename = 'NGC1387_combine_dirty24',
        stokes = 'I',
        outframe = 'LSRK',
        restfreq = '230.538GHz',
        uvrange = '<1750',
        specmode = 'cube',
        nchan = nChanCO,
        start = startCO, 
        width = '2km/s',
        cell=[myCell_dirty],
        imsize=[myImsize_dirty, myImsize_dirty],
        gridder = 'mosaic',
        pbcor = False,
        restoringbeam = '',

        weighting = 'briggs',
        robust = robust,
        # uvtaper=['0.22arcsec'],

        niter = 0,
        chanchunks=-1,
        )

chanstat=imstat(imagename='NGC1387_combine_dirty24.image',chans='0~8')
rms1= chanstat['rms'][0]
chanstat=imstat(imagename='NGC1387_combine_dirty24.image',chans='91~99')
rms2= chanstat['rms'][0]
rms=0.5*(rms1+rms2)
print rms *1000


# clean 1,2,3, test run for velocity range under bary frame (#1 & #2 too large cell size, causing task failure)
# emission from 1208 km/s to 1362 km/s. Choose range 1185 km/s to 1385 km/s
# 4 forgot to change niter from zero
# 5, beam same as dirty13, line-free channels rms=1.10636930437, very similar to dirty13, robust = 1.0, good mask, threshold = '1.5mJy', pbmask = 0.2, sidelobethreshold = 2.0, noisethreshold = 4.25, minbeamfrac =  0.3, lownoisethreshold = 1.2
# 6, robust = 1.5, initial masks too small (actually it's normal), aborted halfway
# 7, robust = 1.5, beam same as dirty10, change lownoisethreshold = 1.0, manual mask for 0-40 channels
# 8, robust = 1.0, niter=10, to construct a mask, forgot to turn down automasking, aborted
# 9, robust=1.0, good mask, niter=10


tclean( vis = 'NGC1387_combine.ms',
        imagename = 'NGC1387_combine_clean9',
        cell=['0.04arcsec'],
        imsize=[750,750], #covers about the FoV 
        gridder = 'mosaic',
        pbcor = True,
        restoringbeam = '',
        specmode = 'cube',
        outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', 
        width = '2km/s', nchan = 100,

        mask='NGC1387_combine_clean9.mask',

        weighting = 'briggs',
        robust = 1.0,
        # uvtaper=uvtaper,

        deconvolver = 'hogbom',
        scales=[0,4,12],# in pixels, equivalent to point source and 1, 3 times the beam.
        niter = 10, # 1000000

        threshold = '1.5mJy',
        interactive = True,
        # usemask = 'auto-multithresh',
        pbmask = 0.2,
        sidelobethreshold = 2.0, # 2.0 recommended
        noisethreshold = 4.25,  # 4.25 recommended
        minbeamfrac =  0.3,       # 0.3 recommended    
        lownoisethreshold = 1.0, # 1.5 recommended
        negativethreshold = 0.0,
        fastnoise=False,

        chanchunks=-1,
        )

# calculate from non-pb corrected, clean cube
chanstat=imstat(imagename='NGC1387_combine_clean5.image',chans='0~8')
rms1= chanstat['rms'][0]
chanstat=imstat(imagename='NGC1387_combine_clean5.image',chans='91~99')
rms2= chanstat['rms'][0]
rms=0.5*(rms1+rms2)
print rms*1000 #


# immoments('NGC1387_combine_clean0.image',
#     moments =[0],
#     axis='spectral',
#     chans = '9~86', # 9~86 for ACA NGC_1387_CO_cube_ACA_clean0.image.pbcor 11~88 for 12-m cube  NGC1387_CO_cube.image.pbcor
#     # excludepix = [-100000.0,0.0],
#     outfile = 'NGC1387_combine_clean0.nonpbcor.moment'
#     )

# immoments('NGC1387_combine_dirty.image',
#     moments =[0],
#     axis='spectral',
#     chans = '9~86', # 9~86 for ACA NGC_1387_CO_cube_ACA_clean0.image.pbcor 11~88 for 12-m cube  NGC1387_CO_cube.image.pbcor
#     # excludepix = [-100000.0,0.0],
#     outfile = 'NGC1387_combine_dirty.nonpbcor.moment'
#     )


# rootfile = 'NGC1387_combine_clean9.mask'
exportfits(
    imagename=rootfile,
    fitsimage=rootfile+'-2.fits',
    dropstokes=True,
    dropdeg=False,
    )

rootfile='/Users/liangf/work/ngc1387/calibrated_all/trial_combine_clean/NGC1387_combine_clean5.image.pbcor'
exportfits(
    imagename=rootfile,
    fitsimage=rootfile+'-VEL.fits',
    dropstokes=True,
    dropdeg=True,
    velocity=True,
    )

# ia.fromfits(outfile = , infile=)

immoments('NGC1387_combine_dirty.image',
    moments =[0],
    axis='spectral',
    chans = '9~86', # 9~86 for ACA NGC_1387_CO_cube_ACA_clean0.image.pbcor 11~88 for 12-m cube  NGC1387_CO_cube.image.pbcor
    # excludepix = [-100000.0,0.0],
    outfile = 'NGC1387_combine_dirty.nonpbcor.moment'
    )

''' Combine 7-m and 12-m data (end) '''



''' Fit the continuum '''

res = imfit(
    imagename = '/Users/liangf/work/continuum-fitting/NGC1387_continuum_2.image.pbcor',
    # excludepix = [-1e10,0],
    residual = '/Users/liangf/work/continuum-fitting/fit-residual',
    model = '/Users/liangf/work/continuum-fitting/fit-model',
    estimates = '/Users/liangf/work/continuum-fitting/fit-guess.txt', # or fit-guess.txt
    logfile = '/Users/liangf/work/continuum-fitting/fit-log.txt',
    newestimates = '/Users/liangf/work/continuum-fitting/fit-guess_new.txt',
    summary = '/Users/liangf/work/continuum-fitting/fit-summary.txt',
    rms = 0.000024,
    )




''' For Jacob's continuum task Elford+2024 '''

myMS = 'uid___A002_Xc44eb5_X2e5d.ms.split.cal.split2'
contSPW_X2e5d_low = '0:0~70;90~127,3:0~1700;2200~3839' # X2e5d
contSPW_X2e5d_up = '1,2'
# 0 232503.990, 1 243962.236, 2 245953.735, 3 230497.653

myMS = 'uid___A002_Xba26cb_X492.ms.split.cal.split2'
contSPW_X492_low = '0,3:0~1700;2200~3839' # X492
contSPW_X492_up = '1:0~20;40~127,2'
# 0 232492.243, 1 243949.806, 2 245941.205, 3 230485.564

myMS = 'uid___A002_Xbbfdf3_X1224.ms.split.cal.split2'
contSPW_X1224_low = '0,3:0~1700;2200~3839' # X1224
contSPW_X1224_up = '1:15~127,2'
# 0 232482.045, 1 243939.017, 2 245930.327, 3 230475.071

myMS = 'uid___A002_Xc2bb44_X2142.ms.split.cal.split2'
contSPW_X2142_low = '0:0~900;1150~2047,1:0~75;90~127'
contSPW_X2142_up = '2,3'
# 0 230580.444, 1 232564.817, 2 244522.322, 3 246514.019


tclean(vis=['uid___A002_Xc44eb5_X2e5d.ms.split.cal.split2',
    'uid___A002_Xba26cb_X492.ms.split.cal.split2',
    'uid___A002_Xbbfdf3_X1224.ms.split.cal.split2',
    'uid___A002_Xc2bb44_X2142.ms.split.cal.split2',],
      imagename = 'NGC1387_upper_continuum.clean',
      spw=[contSPW_X2e5d_up, contSPW_X492_up, contSPW_X1224_up, contSPW_X2142_up],
      threshold='0.09mJy', # measured rms=0.03mJy
      niter=10000, # 10000
      cell=['0.04arcsec'],
      imsize=[750,750],
      outframe='bary',
      restfreq = '230.538GHz',
      deconvolver='multiscale',
      scales=[0,4,12],# inpixels equivalent to point source and 1, 3 times the beam.
      weighting = 'briggs',
      robust=1.0,
      pbcor=True, 
      specmode='mfs',
      restoringbeam='common',
      nterms=1,        
      chanchunks=-1, 
      interactive=True,
      gridder = 'mosaic',
      )


tclean(vis=['uid___A002_Xc44eb5_X2e5d.ms.split.cal.split2',
    'uid___A002_Xba26cb_X492.ms.split.cal.split2',
    'uid___A002_Xbbfdf3_X1224.ms.split.cal.split2',
    'uid___A002_Xc2bb44_X2142.ms.split.cal.split2',],
      imagename = 'NGC1387_lower_continuum.clean2',
      spw=[contSPW_X2e5d_low, contSPW_X492_low, contSPW_X1224_low, contSPW_X2142_low],
      threshold='0.09mJy', # measured rms=0.03mJy
      niter=1000000, # 10000
      cell=['0.04arcsec'],
      imsize=[750,750],
      outframe='bary',
      restfreq = '230.538GHz',
      deconvolver='multiscale',
      scales=[0,4,12],# inpixels equivalent to point source and 1, 3 times the beam.
      weighting = 'briggs',
      robust=1.0,
      pbcor=True, 
      specmode='mfs',
      restoringbeam='common',
      nterms=1,        
      chanchunks=-1, 
      interactive=False,
      mask = 'NGC1387_lower_continuum.clean.mask',
      gridder = 'mosaic',
      )



# only 12-m, low
tclean(vis=['uid___A002_Xc44eb5_X2e5d.ms.split.cal.split2',
    'uid___A002_Xba26cb_X492.ms.split.cal.split2',
    'uid___A002_Xbbfdf3_X1224.ms.split.cal.split2'],
      imagename = 'NGC1387_lower_continuum.clean3',
      spw=[contSPW_X2e5d_low, contSPW_X492_low, contSPW_X1224_low],
      threshold='0.09mJy', # measured rms=0.03mJy
      niter=100000, # 10000
      cell=['0.04arcsec'],
      imsize=[100,100],
      outframe='bary',
      restfreq = '230.538GHz',
      deconvolver='multiscale',
      scales=[0,4,12],# inpixels equivalent to point source and 1, 3 times the beam.
      weighting = 'briggs',
      robust=1.0,
      pbcor=True, 
      specmode='mfs',
      restoringbeam='common',
      nterms=1,        
      # chanchunks=-1, 
      interactive=False,
      mask = 'NGC1387_lower_continuum.clean.mask',
      gridder = 'standard',
      )

# only 12-m, high
tclean(vis=['uid___A002_Xc44eb5_X2e5d.ms.split.cal.split2',
    'uid___A002_Xba26cb_X492.ms.split.cal.split2',
    'uid___A002_Xbbfdf3_X1224.ms.split.cal.split2'],
      imagename = 'NGC1387_upper_continuum.clean3',
      spw=[contSPW_X2e5d_up, contSPW_X492_up, contSPW_X1224_up],
      threshold='0.09mJy', # measured rms=0.03mJy
      niter=100000, # 10000
      cell=['0.04arcsec'],
      imsize=[100,100],
      outframe='bary',
      restfreq = '230.538GHz',
      deconvolver='multiscale',
      scales=[0,4,12],# inpixels equivalent to point source and 1, 3 times the beam.
      weighting = 'briggs',
      robust=1.0,
      pbcor=True, 
      specmode='mfs',
      restoringbeam='common',
      nterms=1,        
      # chanchunks=-1, 
      interactive=False,
      mask = 'NGC1387_lower_continuum.clean.mask',
      gridder = 'standard',
      )

rootfile = 'NGC1387_upper_continuum.clean.image.pbcor'
exportfits(
    imagename=rootfile,
    fitsimage=rootfile+'.fits',
    dropstokes=True,
    dropdeg=False,
    )

''' For Jacob's continuum task (end)'''


''' to investigate flux variation with cleaning depth '''

# RMS = 1.1 mJy/beam

rootfile = '/Users/liangf/work_pub/data/cube_flux_depth/10/NGC1387_combine_clean10m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean9.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 1.0, # uvtaper=uvtaper,
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '1.0mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile+'.image.pbcor.fits', dropstokes=True)

rootfile = '/Users/liangf/work_pub/data/cube_flux_depth/15/NGC1387_combine_clean15m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean9.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 1.0, # uvtaper=uvtaper,
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '1.5mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile+'.image.pbcor.fits', dropstokes=True)

rootfile = '/Users/liangf/work_pub/data/cube_flux_depth/NGC1387_combine_clean20m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean9.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 1.0, # uvtaper=uvtaper,
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '2.0mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile+'.image.pbcor.fits', dropstokes=True)

rootfile = '/Users/liangf/work_pub/data/cube_flux_depth/NGC1387_combine_clean30m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean9.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 1.0, # uvtaper=uvtaper,
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '3.0mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile+'.image.pbcor.fits', dropstokes=True)




''' changing to mask5 '''

rootfile = '/Users/liangf/work_pub/data/cube_flux_depth-mask5/10/NGC1387_combine_clean10m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean5.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 1.0, # uvtaper=uvtaper,
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '1.0mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile+'.image.pbcor.fits', dropstokes=True)
exportfits(imagename=rootfile+'.residual', fitsimage=rootfile+'.residual.fits', dropstokes=True)
exportfits(imagename=rootfile+'.model', fitsimage=rootfile+'.model.fits', dropstokes=True)
# exportfits(imagename=rootfile+'.pb', fitsimage=rootfile+'.pb.fits', dropstokes=True)

# to reproduce the cube in use
rootfile = '/Users/liangf/work_pub/data/cube_flux_depth-mask5/15/NGC1387_combine_clean_paper'
rootfile2 = '/Users/liangf/work_pub/data/cube_flux_depth-mask5/15/NGC1387_combine_clean15m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean5.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 1.0, # uvtaper=uvtaper,
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '1.5mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile2+'.image.pbcor.fits', dropstokes=True)
exportfits(imagename=rootfile+'.residual', fitsimage=rootfile2+'.residual.fits', dropstokes=True)
exportfits(imagename=rootfile+'.model', fitsimage=rootfile2+'.model.fits', dropstokes=True)


rootfile = '/Users/liangf/work_pub/data/cube_flux_depth-mask5/20/NGC1387_combine_clean20m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean5.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 1.0, # uvtaper=uvtaper,
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '2.0mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile+'.image.pbcor.fits', dropstokes=True)
exportfits(imagename=rootfile+'.residual', fitsimage=rootfile+'.residual.fits', dropstokes=True)
exportfits(imagename=rootfile+'.model', fitsimage=rootfile+'.model.fits', dropstokes=True)

rootfile = '/Users/liangf/work_pub/data/cube_flux_depth-mask5/30/NGC1387_combine_clean30m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean5.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 1.0, # uvtaper=uvtaper,
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '3.0mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile+'.image.pbcor.fits', dropstokes=True)
exportfits(imagename=rootfile+'.residual', fitsimage=rootfile+'.residual.fits', dropstokes=True)
exportfits(imagename=rootfile+'.model', fitsimage=rootfile+'.model.fits', dropstokes=True)



''' mask5, low-resolution cube  robust = 2.0, uvtaper='600klambda' '''

rootfile = '/Users/liangf/work_pub/data/cube_flux_depth-mask5-low_res/10/NGC1387_combine_clean10m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean5.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 2.0, uvtaper='600klambda',
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '1.0mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile+'.image.pbcor.fits', dropstokes=True)
exportfits(imagename=rootfile+'.residual', fitsimage=rootfile+'.residual.fits', dropstokes=True)
exportfits(imagename=rootfile+'.model', fitsimage=rootfile+'.model.fits', dropstokes=True)

rootfile = '/Users/liangf/work_pub/data/cube_flux_depth-mask5-low_res/15/NGC1387_combine_clean15m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean5.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 2.0, uvtaper='600klambda',
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '1.5mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile+'.image.pbcor.fits', dropstokes=True)
exportfits(imagename=rootfile+'.residual', fitsimage=rootfile+'.residual.fits', dropstokes=True)
exportfits(imagename=rootfile+'.model', fitsimage=rootfile+'.model.fits', dropstokes=True)


rootfile = '/Users/liangf/work_pub/data/cube_flux_depth-mask5-low_res/20/NGC1387_combine_clean20m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean5.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 2.0, uvtaper='600klambda',
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '2.0mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile+'.image.pbcor.fits', dropstokes=True)
exportfits(imagename=rootfile+'.residual', fitsimage=rootfile+'.residual.fits', dropstokes=True)
exportfits(imagename=rootfile+'.model', fitsimage=rootfile+'.model.fits', dropstokes=True)

rootfile = '/Users/liangf/work_pub/data/cube_flux_depth-mask5-low_res/30/NGC1387_combine_clean30m'
tclean( vis = '/Users/liangf/work_pub/data/NGC1387_combine.ms.contsub',
        mask='/Users/liangf/work_pub/data/NGC1387_combine_clean5.mask', usemask = 'user',
        imagename = rootfile,
        cell=['0.04arcsec'], imsize=[750,750], gridder = 'mosaic', pbcor = True, restoringbeam = '',
        specmode = 'cube', outframe = 'bary', restfreq = '230.538GHz', start = '1185km/s', width = '2km/s', nchan = 100,
        weighting = 'briggs', robust = 2.0, uvtaper='600klambda',
        deconvolver = 'hogbom', niter = 10000000, scales=[0,4,12], # pixels, equivalent to point source / zero and 1, 3 times the beam.
        threshold = '3.0mJy', interactive = False,)
exportfits(imagename=rootfile+'.image.pbcor', fitsimage=rootfile+'.image.pbcor.fits', dropstokes=True)
exportfits(imagename=rootfile+'.residual', fitsimage=rootfile+'.residual.fits', dropstokes=True)
exportfits(imagename=rootfile+'.model', fitsimage=rootfile+'.model.fits', dropstokes=True)


rootfile = '/Users/liangf/work_pub/data/NGC1387_combine_clean5'
exportfits(imagename=rootfile+'.pb', fitsimage=rootfile+'.pb.fits', dropstokes=True)

''' to investigate flux variation with cleaning depth end '''


