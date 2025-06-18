PRO make_galaxy_masktoo, $
	   infile=infile, max_dist = max_dist, $
       ingalpar = ingalpar, outdist = outdist, outmask = outmask

;+
; NAME:
;   make_galaxy_masktoo
;
; PURPOSE:
;   Generate a mask of data cube within distance of max_dist from
;   galaxy center
;
; INPUTS:
;   infile    --   the fits file of original data cube
;   max_dist  --   the maximum distance of mask away from galaxy 
;                  center in [pc]
;   ingalpar     --  the dat file of galaxy parameters
;    
; OUTPUTS:
;   outmask     the fits file of output mask cube
; 
; HISTORY:
;   20190410  LJ   introduced


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

; READ ingalpar
  pa=0.0
  inc = 0.0
  ra_c = ''
  dec_c = ''
  dist_pc = 0.0
  blank = ''
  vsys = 0.0

  openr,lun,ingalpar,/get_lun
  readf,lun,blank
  readf,lun,blank
  readf,lun,pa
  readf,lun,blank
  readf,lun,inc
  readf,lun,blank
  readf,lun,ra_c
  readf,lun,blank
  readf,lun,dec_c
  readf,lun,blank
  readf,lun,dist_pc
  readf,lun,blank
  readf,lun,vsys

  close,lun
  free_lun,lun

  if ((pa eq 0.) and (ra_c eq 0.) and (dec_c eq 0.)) then begin
      print, 'Error: Wrong galaxy name, or No records for the galaxy !'
      return
  endif

; GALAXY CENTER PIXEL
  extast, hdr, astr
  ihr = float((strsplit(ra_c,'h',/extract))[0])
  imin = float((strsplit((strsplit(ra_c,'h',/extract))[1], 'm', /extract))[0])
  xsec = float((strsplit((strsplit(ra_c,'m',/extract))[1], 's', /extract))[0])
  ideg = float((strsplit(dec_c,'d',/extract))[0])
  imn = float((strsplit((strsplit(dec_c,'d',/extract))[1], 'm', /extract))[0])
  xsc = float((strsplit((strsplit(dec_c,'m',/extract))[1], 's', /extract))[0])
  ra_c = ihr * (360./24.) + imin * (360./(24.*60.)) + $
      xsec * (360./(24. * 60. * 60))     ; in [degree]
  dec_c = abs(ideg) + imn * (1./60.) + xsc * (1./3600.)   ; in [degree]
  dec_c = dec_c * ideg/abs(ideg)

  ad2xy, ra_c, dec_c, astr, gal_xc, gal_yc
  ;print,'==========================================='
  print,'Center of Galaxy in Pixel:', gal_xc, gal_yc
  ;print,'==========================================='

; GENERATE MASK
  xy2ad, [0,1], [0,0], astr, ra, dec
  degperpix = sphdist(ra[0], dec[0], ra[1], dec[1], /deg)
  pcperpix= degperpix*!dtor*dist_pc
  axis_ratio = cos(inc/90.*!pi/2.)  ; major and minor axis ratio of gas disk

  mask = cube*0.0
  sz = (size(cube))[1:3]

  dist_cube = mask*0.0 + !values.f_nan

  for i=1, sz[0] do begin
     for j=1,sz[1] do begin
        this_dist = calc_distance(x1=i,y1=j,x2=gal_xc,y2=gal_yc, $
               in_pa = pa, in_inc=inc, pcperpix=pcperpix)  ; in [pc]

		dist_cube[i-1,j-1,*] = this_dist

        if this_dist le max_dist then $
          mask[i-1,j-1,*] = 1
     endfor
  endfor



; WRITE OUT GALAXY 3D MASK FITS FILE
  writefits,outmask,mask,hdr

; WRITE OUT outdist
  if n_elements(outdist) gt 0 then begin
     writefits,outdist,dist_cube,hdr
  endif

END
