#Code to generate v_circ for each LoS:
#Author: MDSmith
#Date: 14/11/17

#Version History: 
# v1: 14/11/17 Constructed for Lijie Liu for NGC4429, MDSmith, Oxford

#This code calculates the line-of-sight circular velocity for each pixel location x,y
#given parameters describing the galaxy mass distribution and geometry

#Note this uses code from Tim Davis' KinMS and Michele Cappellari's JAM packages

# Structure
# main()
# -> calcVCirc_gal()
# -> mge_vcirc(), calculate qintr and dens
# -> vcirc(), quadva integration for integ and then vc

# -> rotation curve table, velocity field fits


from jampy.mge_vcirc import mge_vcirc
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.io import fits
import os
import sys
from astropy.wcs import WCS
from astropy.coordinates import Angle

def plotMap(xPos,yPos,quant,title = '',show=True):
    ###############################################################
    #Place interesting quantity into array:
    import matplotlib.pyplot as plt
    print('max(xPos),max(yPos)',max(xPos),max(yPos))
    xsize = int(max(xPos))
    ysize = int(max(yPos))
    outputArray = np.zeros([xsize, ysize])
    cellIndex = np.arange(0,xsize-1)
    cent=np.array([0,0])
    for i in np.arange(0,xPos.shape[0]):
        cubeX=int(xPos[i]+cent[0])
        cubeY=int(yPos[i]+cent[1])
        #print (cubeX,cubeY)
        if (0 < cubeX < (xsize-1)) & (0 < cubeY < (ysize-1)):
            #print 'Adding cloud to output cube'
            outputArray[cubeX,cubeY]=quant[i]
    #plt.contour(cellIndex,cellIndex,outputArray.T,levels=levels)
    # outputArray.T[outputArray.T<-100] = 0 # to account for the central pixel dominated by the SMBH

    plt.figure()
    plt.imshow(outputArray.T,cmap='jet',origin='lower',vmin =-100, vmax=100) # edited by Eric
    # print('max',outputArray[int(xsize/2-10):int(xsize/2 +10),int(xsize/2-10):int(xsize/2 +10)])
    plt.title(title)
    plt.colorbar()
    plt.savefig('/Users/ericliang/n1387/work_pub/plot/circular_velocity/overlaySSvLos.pdf',transparent=True)
    if show==True:
        plt.show()
    plt.close()
    ################################################################


def calcRadius(xPos,yPos, posAng, inc,cen =[0.,0.]):
    #Given deprojected distances of a set of positions [xPos,yPos] from cen
    #in unit of xpos/pos/cen 
    ang2Rot = 90 - posAng      # in [degree]
    c = np.cos(np.radians(ang2Rot))
    s = np.sin(np.radians(ang2Rot))
    delX=xPos-np.full(xPos.shape,cen[0])
    delY=yPos-np.full(yPos.shape,cen[1])
    majComp = delX * c - delY * s
    minComp = delX * s + delY * c
    #De-project using the inclination
    minProj = minComp/np.cos(np.radians(inc))
    
    r_flat = np.sqrt((majComp * majComp) + (minProj * minProj))
    return r_flat

def kinms_create_velField_oneSided(velRad,velProf,r_flat,inc,posAng,xPos,yPos,vPhaseCent=[0.0,0.0],vPosAng=False,vRadial=0.0,posAng_rad=0.0,inc_rad=0.0, sbMode = "axisymmetric",cloudIntensity=None):
    """
    Function modified from Tim Davis' KinMSpy by MDSmith
    
    This function takes the input circular velocity distribution
    and the position of point sources and creates the velocity field 
    taking into account warps, inflow/outflow etc as required.
    
    Parameters
    ----------
    velRad : np.ndarray of double
            Radius vector (in units of pixels or arcsec).
    
    velProf : np.ndarray of double
            Velocity profile (in units of km/s) at velRad.
    
    r_flat : np.ndarray of double
            Radius of each cloudlet from the kinematic centre
            in the plane of the disc. Units of pixels or arcsec (same unit
            as velRad).
            
    inc : double or np.ndarray of double in [degree]
            Inclination of the disc, using the usual astronomical convention.
            Can be either a double, or an array of doubles. If single valued
            then the disc is flat. If an array is passed then it should
            describe how the galaxy inclination changes as a function of `velrad`.
            Used to create inclination warps.
        
    posAng : double or np.ndarray of double in [degree]
            Position angle of the disc, using the usual astronomical convention.
            Can be either a double, or an array of doubles. If single valued
            then the disc major axis is straight. If an array is passed then it should
            describe how the position angle changes as a function of `velrad`.
            Used to create position angle warps.
        
    gasSigma : double or np.ndarray of double in [km/s]
            Velocity dispersion of the gas. Units of km/s. 
            Can be either a double, or an array of doubles. If single valued
            then the velocity dispersion is constant throughout the disc.
            If an array is passed then it should describe how the velocity
            dispersion changes as a function of `velrad`.
        
    seed : list of int
            List of length 4 containing the seeds for random number generation.
    
    xPos : np.ndarray of double in [pixel] in final output coordinate
            X position of each cloudlet. Units of pixels. 
    
    yPos : np.ndarray of double in [pixel] in final output coordinate
            Y position of each cloudlet. Units of pixels. 
    
    vPhaseCent : list of double in [pixel, pixel]
         (Default value = [0, 0])
            Kinematic centre of the rotation in the x-y plane. Units of pixels.
            Used if the kinematic and morphological centres are not the same.
    
    vPosAng : double or np.ndarray of double  in [degree]
         (Default value = False)
            Kinematic position angle of the disc, using the usual astronomical convention.
            Can be either a double, or an array of doubles. If single valued
            then the disc kinematic major axis is straight. If an array is passed then it should
            describe how the kinematic position angle changes as a function of `velrad`.
            Used if the kinematic and morphological position angles are not the same.
    
    vRadial : double or np.ndarray of double
         (Default value = 0)
            Magnitude of inflow/outflowing motions (km/s). Negative
            numbers here are inflow, positive numbers denote
            outflow. These are included in the velocity field using
            formalism of KINEMETRY (Krajnovic et al. 2006 MNRAS, 366, 787). 
            Can input a constant or a vector, giving the radial
            motion as a function of the radius vector
            `velrad`. Default is no inflow/outflow.
    
    posAng_rad : double or np.ndarray of double
         (Default value = 0)
            Position angle of the disc at the position `r_flat` of each cloudlet.
    
    inc_rad : double or np.ndarray of double
         (Default value = 0)
            Inclination angle of the disc at the position `r_flat` of each cloudlet.
    
    Returns
    -------
    los_vel : np.ndarray of double
            Line of sight velocity of each cloudlet, in km/s.
    
    """
    inc_rad = np.full(xPos.shape,inc)
    velInterFunc = interpolate.interp1d(velRad,velProf,kind='linear')
    #print np.max(r_flat) 
    #print np.min(r_flat)
    # vRad - linear rotation velocity at each pixel of map in [km/s]
    vRad = velInterFunc(r_flat)
    # los_vel - line of sight velocity at each pixel of map in [km/s]
    los_vel = np.empty(len(vRad))
    
    # Find the rotation angle so the velocity field has the correct kinematic position angle (allows warps)
    if sbMode == 'axisymmetric':
        if not vPosAng:
            ang2rot=0.0
        else:
            if isinstance(vPosAng, (list, tuple, np.ndarray)):
                vPosAngInterFunc = interpolate.interp1d(velRad,vPosAng,kind='linear')
                vPosAng_rad = vPosAngInterFunc(r_flat)
            else:
                vPosAng_rad = np.full(len(r_flat),vPosAng,np.double)
            ang2rot = ((posAng_rad-vPosAng_rad))
    elif sbMode=='skyProfile':
        #print 'Calculating los velocities in skySampler mode'
        #Includes rotation to the correct kinematic major axis, as this cannot be performed by the global rotation as in the normal KinMS code
        #vPosAng_rad: position angle of kinematic axis at each r_flat in
        #[degree]
        if not vPosAng:
            vPosAng_rad = np.full(len(r_flat),0.,np.double)
        else:
            if isinstance(vPosAng, (list, tuple, np.ndarray)):
                vPosAngInterFunc = interpolate.interp1d(velRad,vPosAng,kind='linear')
                vPosAng_rad = vPosAngInterFunc(r_flat)
            else:
                vPosAng_rad = np.full(len(r_flat),vPosAng,np.double)
            ang2rot = ((90. + posAng_rad-vPosAng_rad))
    #Calculate the correct tan components - erroneous 'W' plots are due to incorrect combination of coTan and vRad, which produces the velocity field
    
    projFactor = np.sin(np.radians(inc_rad))
    #print inc_rad
    if sbMode == 'axisymmetric':
        coTan = np.cos(np.arctan2((yPos),(xPos))+ (np.radians(ang2rot)))
    elif sbMode == 'skyProfile':
        vPosAng = posAng
        # rotation angle (counterclock) of new coordinate from old coordinate
        kPArot = vPosAng_rad - 90   
        cPA = np.cos(np.radians(kPArot))
        sPA = np.sin(np.radians(kPArot))
        delX = xPos - vPhaseCent[0]
        delY = yPos - vPhaseCent[1]
        #majComp = -sPA * delX + cPA *delY
        #minComp = cPA * delX + sPA * delY
        #c=np.full(minComp.shape,np.cos(np.radians(inc)))    
        #minProj = np.multiply(minComp,c)
        #coTan = np.cos(np.arctan2((majComp),(minProj)))

        majComp = delX*cPA + delY*sPA   # deprojected maj component of [xpos, ypos]
        minComp = -delX*sPA + delY*cPA   
        c=np.full(minComp.shape,np.cos(np.radians(inc)))    
        minComp = minComp/c            # deprojected min component of [xpos, ypos]
        
        coTan = np.cos(np.arctan2((minComp),(majComp)))
    #Calculate the los velocity for each cloudlet
    los_vel = (-1) * vRad * coTan * projFactor           
    
    
    #Add radial inflow/outflow
    if vRadial != 0:
        if isinstance(vRadial, (list, tuple, np.ndarray)):
            vRadialInterFunc = interpolate.interp1d(velRad,vRadial,kind='linear')
            vRadial_rad = vRadialInterFunc(r_flat)
        else:
            vRadial_rad=np.full(len(r_flat),vRadial,np.double)
        los_vel += vRadial_rad * coTan * projFactor
    
    return los_vel, coTan, vRad, projFactor


def makePsi(velrad, psi, debug=False):
    # Returns the appropriate M/L vector based on input
    psi = np.zeros(velrad.shape)
    
    # Case 1: constant M/L ratio
    if len(ml_par) == 1:
        psi = psi + ml_par[0]
    
    if len(ml_par) == 4:
        psiCent = ml_par[0]
        r_break = ml_par[1]
        psiOuter = ml_par[2]
        rOuter = ml_par[3]

        for i in np.arange(0, velrad.shape[0]):
            r = velrad[i]
            if r <= r_break:
                psi[i] = psiCent
            elif (r_break < r < rOuter):
                m2 = (psiOuter - psiCent) / (rOuter - r_break)
                c2 = psiCent - m2 * r_break
                psi[i] = m2 * r + c2
            else:
                psi[i] = psiOuter
    
    if debug:
        # plot the generated function for the user
        plt.plot(velrad, psi)
        plt.xlim(0., 1.5 * rOuter)
        plt.ylim(0., psiOuter * 3.)
        plt.xlabel('Radius')
        plt.ylabel('M/L')
        plt.show()
        plt.close()

    return psi
    
def calcLoSVel():
    posAng = 90 - posAng    #Corrects for PA oriented aclockwise from celestial north (the y axis, but points to galaxy x axis)
    if isinstance(posAng, (list, tuple, np.ndarray)):
        posAngRadInterFunc = interpolate.interp1d(velRad,posAng,kind='linear')
        posAng_rad = posAngRadInterFunc(r_flat)
    else:
        posAng_rad = np.full(len(r_flat),posAng,np.double)
    
    if isinstance(inc, (list, tuple, np.ndarray)):
        incRadInterFunc = interpolate.interp1d(velRad,inc,kind='linear')
        inc_rad = incRadInterFunc(r_flat)
    else:
        inc_rad = np.full(len(r_flat),inc,np.double)


def calcVCirc_gal(datacube='', distance=0., inc=0., posAng=0.,\
        cen = 0., vsys = 0., mge_file='', ml_par = [], MBH = 0., \
        rotdat = '', rotfits='', setunit = ''):
    # 0. Clean up
    if os.path.exists(rotdat):
         os.remove(rotdat)
         print(rotdat+' exists')
         # sys.exit()
    if os.path.exists(rotfits):
         os.remove(rotfits)
         print(rotfits+' exists')
         # sys.exit()
    if os.path.exists('overlaySSvLos.pdf'):
         os.remove('overlaySSvLos.pdf')
         print('overlaySSvLos.pdf exists')
         # sys.exit()

    
    # 2. Read data cube header
    dat = fits.open(datacube)
    linecube = dat[0].data
    hdr = dat[0].header
    cellsize = abs(hdr['CDELT1'])             # in [degree]
    cellsize = cellsize * 3600.               # in [arcsec]
    BMAJ = hdr['BMAJ']
    BMIN = hdr['BMIN'] 
    BPA = hdr['BPA']
    beam = np.array([BMAJ*3600., BMIN*3600., BPA])  # in [arcsec]
    NAXIS1 = hdr['NAXIS1']
    NAXIS2 = hdr['NAXIS2']
    NAXIS3 = hdr['NAXIS3']
    crval3 = hdr['CRVAL3']
    cdelt3 = hdr['CDELT3']
    crpix3 = hdr['CRPIX3']
    CUNIT3 = hdr['CUNIT3']
    mapSize = [NAXIS1, NAXIS2]    

    
    # 3. MGC parameterisation of galaxy light profile from Davis et al (2017)
    data = np.loadtxt(mge_file, delimiter=",", comments='#')
    #   V-band surface brightness in L_solar pc^-2  # potentially modified by Eric
    mge_comp = data[:, 0]
    #   Stardard deviation (width) in [arcsec]
    mge_width = data[:, 1]
    #   Axis ratio
    mge_q = data[:, 2]
    
    # 4. Velocity field by stellar mass
    #   Setup radial vector to discretise v_circ onto
    rad = np.arange(0., abs(NAXIS1*cellsize/np.cos(np.radians(inc))),0.001)   # in [arcsec]
    rad[0] = 1e-8
    
    # print ( [mge_comp, mge_width, mge_q, inc, 0., distance, rad] )
    # print ( [np.array([0.]), np.array([1.]), np.array([1.]), inc, MBH, distance, rad] )
    # sys.exit()

    #   4.1 Assume M/L = 1 since this galaxy has variable M/L
    #   mge_vcirc(surf_pc, sigma_arcsec, qobs, inc_deg, mbh, distance, rad, soft=0.)
    velProf = mge_vcirc(mge_comp, mge_width, mge_q, inc, 0., distance, rad)    

    # # optional plot
    # plt.figure()
    # plt.plot(rad,velProf)
    # plt.gca().set_ylim([0.,max(velProf)*1.3])
    # plt.gca().set_xlim([0.,max(rad)*1.3])
    # plt.xlabel('R [arcsec]')
    # plt.ylabel('v [km/s]')
    # plt.title('Rotation curve assuming ml = 1')
    # plt.show()
    # plt.close()
    
    #   4.2 Scales the velocity profile to the variable mass-to-light ratio 
    #   makePsi(velrad, psiCent, r_break_arcsec, psiOuter, rOuter_arcsec, debug=False):
    psi = makePsi(rad, ml_par)
    velProf = np.multiply(velProf,np.sqrt(psi))  
    velProf_star = velProf.copy()
    
    # 5. Velocity field by black hole and total velocity field
    velBH = mge_vcirc(np.array([0.]), np.array([1.]), np.array([1.]), inc, MBH, distance, rad)
    
    velProf = np.sqrt(velBH ** 2 + velProf ** 2)

    fid = open(rotdat, 'w')
    rad_pc = rad/(180.*3600.)*np.pi*distance*1e6  # [pc]
    fid.write(';rad_pc, circular_velocity_kms\n')
    thedata = np.column_stack((rad_pc, velProf))
    np.savetxt(fid, thedata, fmt=['%1.4e']*2,delimiter=",")
    fid.close()

    # 6. Generate rotation fits file
    #   6.1 Make array of coordinates for each pixel
    xPos = np.zeros(mapSize[0]*mapSize[1])
    yPos = np.zeros(mapSize[0]*mapSize[1])
    k=0
    for i in np.arange(0,mapSize[0]):
        for j in np.arange(0,mapSize[1]):
            xPos[k] = i
            yPos[k] = j
            k=k+1
    
    #      deprojected distance from xPos/yPos to cen in [pixel]
    r_flat = calcRadius(xPos,yPos,posAng,inc,cen)   # in [pixel]
    r_flat = r_flat * cellsize  # in [arcsec]
    
    
    #   6.2 Calculate the losVel for each location (x,y) 
    #       Note we neglect the effect of gas mass
    rad[0] = 0.
    los_vel,coTan,vRad,projFactor = kinms_create_velField_oneSided(rad,velProf,r_flat,inc,0.,xPos,yPos,vPosAng = posAng, vPhaseCent = cen,sbMode = "skyProfile",cloudIntensity=None)

    VelArray = np.zeros(mapSize)
    for k in np.arange(0,xPos.shape[0]):
        VelArray[int(xPos[k]),int(yPos[k])] = los_vel[k]
    
    plotMap(xPos,yPos,los_vel,title = 'Velocity Field',show=False) # optional plot

    #   Write VelArray to fits
    rotVfield = np.transpose(VelArray)

    VelCube = np.repeat(rotVfield[np.newaxis,:,:], NAXIS3, axis=0)

    if setunit in ['km/s', 'KM/S']:
        VelCube = VelCube + vsys

        hdu = fits.PrimaryHDU()
        hdu.data = VelCube
        hdu.header
        hdu.header = hdr
        hdu.header['BTYPE'] = 'Velocity'
        hdu.header['BUNIT'] = 'km/s'
        hdu.header['NAXIS1'] = NAXIS2
        hdu.header['NAXIS2'] = NAXIS1
    
        hdu.writeto(rotfits)
    elif setunit in ['pixel','channel']:
        if (CUNIT3.replace(' ','')) in ['m/s','M/S']:
            crval3 = crval3/1e3
            cdelt3 = cdelt3/1e3
        #print 'crval3, cdelt3, crpix3',crval3, cdelt3, crpix3
        VelCube = (VelCube +vsys - (crval3))/cdelt3 + crpix3 -1  # in pixel or channel
        hdu = fits.PrimaryHDU()
        hdu.data = VelCube
        hdu.header
        hdu.header = hdr
        #hdu.header.remove('CRVAL3')
        #hdu.header.remove('CDELT3')
        #hdu.header.remove('CUNIT3')
        #hdu.header.remove('CTYPE3')
        hdu.header['BTYPE'] = 'Velocity'
        hdu.header['BUNIT'] = 'CHANNEL'
        hdu.header['NAXIS1'] = NAXIS2
        hdu.header['NAXIS2'] = NAXIS1
    
        hdu.writeto(rotfits)


    # 8. Plot rotation curves
    plt.figure()
    plt.plot(rad,velProf,label='Total')
    plt.plot(rad,velProf_star,ls='--',label=r'Stellar ($\Upsilon$ = %.2f)'%ml_par[0])
    plt.plot(rad,velBH,ls='--',label='SMBH')
    plt.legend()
    plt.gca().set_ylim([0.,400.])
    plt.gca().set_xlim([0.,10.])
    plt.title('Circular Velocity Curve')
    plt.xlabel('R [arcsec]')
    plt.ylabel('v [km/s]')
    plt.savefig('/Users/ericliang/n1387/work_pub/plot/circular_velocity/NGC1387_velocity_curve.png')
    # plt.show()
    plt.close()

    rad[0]= 1e-8
    plt.figure()
    plt.plot(rad,velProf/rad)
    plt.gca().set_ylim([0.,400.])
    plt.gca().set_xlim([0.,10.])
    plt.xlabel('R [arcsec]')
    plt.ylabel(r'$\omega$ [km/s/arcsec]')
    plt.title('Angular Circular Velocity Curve')
    plt.savefig('/Users/ericliang/n1387/work_pub/plot/circular_velocity/NGC1387_angular_velocity_curve.png')
    # plt.show()
    plt.close()
    rad[0] = 0.
    
    A=rad * 0.
    for i in np.arange(0,rad.shape[0]-1):
        A[i] = 0.5 * (velProf[i]/rad[i] - (velProf[i+1]-velProf[i])/(rad[i+1]-rad[i]))
    
    plt.plot(rad,A)
    plt.xlabel('R [arcsec]')
    plt.ylabel('A [km/s/\"]')
    plt.ylim(0.,100.)
    plt.savefig('/Users/ericliang/n1387/work_pub/plot/circular_velocity/NGC1387_shear.png')
    # plt.show()
    plt.close()

    A=rad * 0.
    for i in np.arange(0,rad.shape[0]-1):
        A[i] = 0.5 * (velProf[i]/rad_pc[i] - (velProf[i+1]-velProf[i])/(rad_pc[i+1]-rad_pc[i]))
    
    plt.plot(rad_pc,A)
    plt.xlabel('R [pc]')
    plt.ylabel('A [km/s/pc]')
    # plt.ylim(0.,100.)
    plt.xlim(0,1000)
    plt.ylim(0.1,300)
    plt.yscale('log')
    plt.savefig('/Users/ericliang/n1387/work_pub/plot/circular_velocity/NGC1387_shear-pc.png')
    # plt.show()
    plt.close()
    print(rad_pc)
    print(A)

if __name__=="__main__":
    """
    """
    # ========= Only Need to Change Below ========================
    # Input
    gal_property = np.genfromtxt('/Users/ericliang/n1387/work_pub/data/NGC1387_galpar.dat', comments=';',dtype=str)
    posAng = float(gal_property[0])     # in [degree]
    inc = float(gal_property[1])        # in [degree]
    ra = gal_property[2]
    dec = gal_property[3]
    distance = float(gal_property[4]) / 1e6      # in [Mpc]
    vsys = float(gal_property[5])        # in [km/s]

    
    datacube = '/Users/ericliang/n1387/work_pub/data/NGC1387_combine_clean5.image.pbcor.fits'
    hdul=fits.open(datacube)
    ml_par = [float(gal_property[7])]        # mass to light ratio, an array
    MBH = np.power(10, float(gal_property[8]))       # M_sol
    mge_file = '/Users/ericliang/n1387/work_pub/data/NGC1387_mge_gauss.csv'

    wcs = WCS(header=hdul[0].header)
    cen_raw = wcs.all_world2pix([Angle(ra).value*15],[Angle(dec).value],[229], 0, ra_dec_order=True) #
    # cen = [381, 378]    # pixel indices RA = 380.589, DEC=378.1115
    cen = [cen_raw[0][0], cen_raw[1][0]]

    # Output
    # 1.  rotation curve data file
    rotdat = '/Users/ericliang/n1387/work_pub/data/NGC1387_velocity_curve-2.dat'
    # 2.  rotation fits same as mom1 map but it's 3D fits, 
    #     every dimension is same
    rotfits = '/Users/ericliang/n1387/work_pub/data/NGC1387_rotation_vfield_2kms-2.fits'
    # ============================================================
    
    calcVCirc_gal(datacube=datacube,\
        distance=distance, inc=inc, posAng = posAng, cen = cen,\
        vsys=vsys, mge_file=mge_file, ml_par=ml_par, MBH=MBH,\
        rotdat=rotdat, rotfits=rotfits,setunit='km/s')

