FUNCTION PADDING,im,edge

; if edge>0 pad data with some zero-value pixels in each dimension 
; if edge<0 remove the padding pixels

if  edge gt 0 then begin
    struct=size(im, /STRUCTURE)
    imdims=size(im)
    tmp=MAKE_ARRAY(struct.dimensions[0:imdims[0]-1]+2*edge, TYPE=struct.TYPE)
    if imdims[0] eq 1 then tmp[edge]=im
    if imdims[0] eq 2 then tmp[edge,edge]=im
    if imdims[0] eq 3 then tmp[edge,edge,edge]=im
    if imdims[0] eq 4 then tmp[edge,edge,edge,edge]=im
endif
if  edge lt 0 then begin
    imdims=size(im)
    if imdims[0] eq 1 then tmp=im[-edge:(edge+imdims[1]-1)]
    if imdims[0] eq 2 then tmp=im[-edge:(edge+imdims[1]-1),-edge:(edge+imdims[2]-1)]
    if imdims[0] eq 3 then tmp=im[-edge:(edge+imdims[1]-1),-edge:(edge+imdims[2]-1),-edge:(edge+imdims[3]-1)]
    if imdims[0] eq 4 then tmp=im[-edge:(edge+imdims[1]-1),-edge:(edge+imdims[2]-1),-edge:(edge+imdims[3]-1),-edge:(edge+imdims[4]-1)]
endif
if  edge eq 0 then tmp=im

return,tmp
END


PRO make_galaxy_mask, infile=infile, rmsfile=rmsfile,$
	            outmask=outmask, outfile=outfile, $
	            fbeam=fbeam, svel=svel , $
                hi_thresh=hi_thresh, lo_thresh=lo_thresh, $
                chmin=chmin, guard=guard

;+
; NAME:
;   make_galaxy_mask
;
; PURPOSE:
;   Smooth FITs cube/image to lower spatial and/or velosity resolution
;   to generate a mask for the region of galaxy only, and then [option]
;   extract a subimage from the original (before smooth) data cube/image
;   without considering the useless edges and update astrometry 
;   in FITS header
;
; INPUTS:
;   infile    --   the fits file of original data cube
;   rmsfile   --   the fits file of original rms cube (same size as infile)
;   outmask   --   the fits file of output mask cube
;   [outfile] --   if set, then extract a subimage as a output fits file
;   [fbeam]   --   desired beam(psf) size in arcsec (could be a scalar or 3-element vector)
;             --   e.g.  fbeam=5             5"X5" (0D)    2D Gaussian
;             --         fbeam=[10,5,90]     10"x5"(+90d)  2D Gaussian
;   [svel]    --   (for 3d case only) desired spectra beam in km/s, smooth spectra
;                  with a Gaussian function, if svel=0.0, no spectra 
;                  smoothing is performed 
;   hi_thresh --   the minmum snr that required by at least one pixel
;                  in smoothed_cube
;   lo_thresh --   the minmum snr that required by all pixels
;                  in smoothed_cube
;   [chmin]   --   minimum number of velocity channels in initial mask, 
;                  default value is 0
;   [guard]   --   interger describing the number of pixels along each
;                  direction that the mask should be expanded.
;                  A vector integers (guard=[xshift,yshift,vshift]) describing 
;                  the number of pixels along each axis that the mask
;                  should be expaned. If guard has only one element, then
;                  it is for v/z direction, if it has two elements then it 
;                  is for x/y directions.
;    
; OUTPUTS:
;   outmask     the fits file of output mask cube
;   [outfile]   if set, then extract a subimage as a output fits file
;   epsfile     1) mom0_before_and_after_mask.eps;
;               2) smooth_before_and_after.eps
; 
; HISTORY:
;   20161205  LJ   introduced

; DEFAULT INPUT PARAMETERS
if not keyword_set(guard) then guard=0

; READ IN THE DATA AND NOISE CUBE
  if n_elements(infile) gt 0 then begin
     file_data = file_search(infile, count=file_ct)
     if file_ct eq 0 then begin
        message, "Data not found.", /info
        return
     endif else begin
        cube = readfits(file_data, hdr)
     endelse
  endif

  if n_elements(rmsfile) gt 0 then begin
     file_rms = file_search(rmsfile, count=file_ct)
     if file_ct eq 0 then begin
        message, "Noise not found.", /info
        return
     endif else begin
        rms = readfits(file_rms, rms_hdr)
     endelse
  endif

; GET THE ORIGINAL BEAM SIZE
  ipsf=[0.,0.,0.]
  RADIOHEAD,hdr,s=h0
  bunit = SXPAR(hdr, 'BUNIT')
  ipsf[0]=h0.bmaj  ; in arcsec
  ipsf[1]=h0.bmin  ; in arcsec
  ipsf[2]=h0.bpa   ; in degrees (astro convention)

; SET THE DEFAULT FBEAM SIZE TO 3 TIMES LARGER OF ORIGINAL BEAM SIZE
  if n_elements(fbeam) eq 0 then begin
     fbeam=[ipsf[0]*3.,ipsf[1]*3.,ipsf[2]]
  endif

; SMOOTH THE CUBE 
  smooth3d, cube/rms, hdr, im_smo,hd_smo,fbeam,svel=svel,ifail=ifail
  if ifail ne 0 then begin
      errmsg = 'ERROR - the target beam is too small!'
      message,'ERROR - ' + errmsg,level=1
  endif

; CALCULATE SMOOTHED RMS AFTER CLIPPING BASED ON MIRIAD SIGEST.FOR
  nchan=(size(cube,/d))[2]
  sen=err_cube(cube/rms,planes=indgen(nchan))
  mask1 = abs(im_smo) lt 2.5*sqrt(!PI/2)*meanabsdev(im_smo,/nan)
  sen_smo=err_cube(im_smo,planes=indgen(nchan),mask=mask1)

; GENERATE 2D IMAGES AND FITS FOR CUBE/RMS AND IM_SMO
  im1=total(cube/rms,3,/nan)
  SXADDPAR,hdr,'DATAMAX', max(im1,/nan), before='HISTORY'
  SXADDPAR,hdr,'DATAMIN', min(im1,/nan), before='HISTORY'
  SXADDPAR,hdr,'BUNIT', 'sigma', before='HISTORY'
  writefits,'sigma_before_smooth.fits',im1,hdr
  
  im2=total(im_smo,3,/nan)
  SXADDPAR,hd_smo,'DATAMAX', max(im2,/nan), before='HISTORY'
  SXADDPAR,hd_smo,'DATAMIN', min(im2,/nan), before='HISTORY'
  SXADDPAR,hd_smo,'BUNIT', 'sigma', before='HISTORY'
  writefits,'sigma_after_smooth.fits',im2,hd_smo
  
  ;plt2d, infile1='sigma_before_smooth.fits', $
  ;	     infile2='sigma_after_smooth.fits',$
  ;		 outfile='smooth_before_and_after.eps',$
  ;	     label={title1:'Intensity', title2:'Intensity'}

  file_delete,'sigma_before_smooth.fits'
  file_delete,'sigma_after_smooth.fits'

; GENERATE THE 3D MASK FOR WHOLE IMAGE THAT SATISFY HI_THRESH
; LO_THRESH AND CHMIN                                         
  if  hi_thresh gt 0.0 then begin
      mask = cube*0.0
      mask[where(im_smo gt sen_smo*hi_thresh, ct, /null)]=1.0
      ; REQUIRE MINIMUM NUMBER OF CHANNELS
      if chmin gt 1 then begin
          mask = padding(mask,chmin-1)
          for j = 2, chmin do begin
              shmask = mask*0.0
              for k = 1, j-1 do begin
                  shmask += shift(mask,0,0,k) + shift(mask,0,0,-k)
              endfor
              mask = mask * (shmask-j+2) gt 0
          endfor
      endif
      ; DILATE THE MASK IF REQUESTED
      if (total(mask,/nan) gt 0.0) and (lo_thresh gt 0.0) then begin
          ; CONSTRAINT MASK
          lo_threshmask = im_smo gt lo_thresh*sen_smo
          ; REQUIRE MINIMUM NUMBER OF CHANNELS (CURRENTLY DISABLED)
          if chmin gt 1 then begin
              lo_threshmask = padding(lo_threshmask,chmin-1)
          endif
          ; EXPAND MASK
          mask = float(dilate_mask(mask, constraint = lo_threshmask))
      endif
      if chmin gt 1 then mask = padding(mask,-chmin+1)
      ; ADD A GUARD BAND IF REQUESTED
      if total(guard) gt 0 then begin
          mask = maskguard_3d(mask, guard = guard)
    endif
  endif else begin
      mask = cube*0.0 + 1.0
  endelse

; MASK ONLY THE REGION FOR GALAXY
  regions = label_region(mask gt 0, all_neighbors=all_neighbors, /ulong)
  if total(regions) eq 0 then begin
    print,'**********************'
	print, 'label region failed'
    print,'**********************'
	return
  endif	
  h = histogram(regions, binsize=1, min=1, rev=ri)
  id=where(h eq max(h))
  inds=ri[ri[id]:(ri[id+1])-1]
  galaxy_mask=cube*0.0
  galaxy_mask(inds) = 1
 
  ; VECTORIZE inds
  ind_to_xyv, inds, x=x, y=y, v=v, sz=size(cube)

; APPLY GALAXY_MASK TO ORIGINAL DATA CUBE
  ;cube2=cube*0.0+!values.f_nan
  cube2=cube*0.0
  cube2[inds]=cube[inds]
  im1=total(cube,3,/nan)*abs(h0.cdelt[2])
  im2=total(cube2,3,/nan)*abs(h0.cdelt[2])

  SXADDPAR,hdr,'DATAMAX', max(im1,/nan), before='HISTORY'
  SXADDPAR,hdr,'DATAMIN', min(im1,/nan), before='HISTORY'
  SXADDPAR,hdr,'BUNIT',strtrim(bunit,2)+'.KM/S', before='HISTORY'
  writefits,'mom0_before_mask.fits',im1,hdr
  
  SXADDPAR,hd_smo,'DATAMAX', max(im2,/nan), before='HISTORY'
  SXADDPAR,hd_smo,'DATAMIN', min(im2,/nan), before='HISTORY'
  SXADDPAR,hd_smo,'BUNIT',strtrim(bunit,2)+'.KM/S', before='HISTORY'
  writefits,'mom0_after_mask.fits',im2,hd_smo
  
  ;plt2d, infile1='mom0_before_mask.fits', $
  ;	     infile2='mom0_after_mask.fits',$
  ;		 outfile='mom0_before_and_after_mask.eps',$
  ;	     label={title1:'Intensity', title2:'Intensity'}

  file_delete,'mom0_before_mask.fits'
  file_delete,'mom0_after_mask.fits'

; WRITE OUT GALAXY 3D MASK FITS FILE
  writefits,outmask,galaxy_mask,hdr

END
