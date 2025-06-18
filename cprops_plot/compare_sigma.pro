PRO compare_sigma, infile = infile $
						 , ingalpar = ingalpar $
						 , inprop = inprop $
					     , rotdat = rotdat $
	                     , outtable=outtable

						 ; , rotStddat = rotStddat $
						 ; , pex_table = pex_table  $
	           ;           , out_vrms_compare=out_vrms_compare $
						 ; , inner_pc = inner_pc $
             ;             , ring_pc = ring_pc $
             ;             , outer_pc = outer_pc $
						 ; , veff_mode = veff_mode $
	           ;           , out_beta_model_fig = out_beta_model_fig $
	           ;           , out_beta_obser_fig = out_beta_obser_fig $
						 ; , logbeta = logbeta $
					   ;   , out_shear_r = out_shear_r


; READ INPUT TABLE
  ; readcol,intable, ID, RA, Dec, Vlsr, R, d_R, vrms, d_vrms, vturb, d_vturb, $
  ;     Lum, d_Lum, Mlum, d_Mlum, Tmax, Omega_shear, Theta_shear, distance, $
  ;     F='I, a, a, f, f, f, f, f, f, f, f, f, f, f, f, f, f, f'
  ; Lum = Lum * 1e4
  ; d_Lum = d_Lum * 1e4
  ; Mlum = Mlum * 1e5
  ; d_Mlum = d_Mlum * 1e5

  restore, inprop

  distance = props.r_gal
  R = props.RADRMS_EXTRAP_DECONV
  D_R = props.RADRMS_EXTRAP_DECONV_UC * R
  ID = props.peaknum

  nclumps = n_elements(distance)

; READ INGALPAR
  if n_elements(ingalpar) gt 0 then begin
; Read galaxy parameters from ingalpar file
  pa=0.0
  inc = 0.0
  ra_c = ''
  dec_c = ''
  dist_pc = 0.0
  blank = ''
  vsys = 0.0
  inner=300.
  outer=600.
  alpha = 0.

  openr,lun,ingalpar,/get_lun
  readf,lun,blank
  readf,lun,blank
  readf,lun,blank
  readf,lun,blank
  readf,lun,pa
  readf,lun,blank
  readf,lun,inc
  readf,lun,blank
  readf,lun,blank
  readf,lun,blank
  readf,lun,ra_c
  readf,lun,blank
  readf,lun,dec_c
  readf,lun,blank
  readf,lun,dist_pc
  readf,lun,blank
  readf,lun,vsys
  readf,lun,blank
  readf,lun,alpha

  close,lun
  free_lun,lun
	
	
	  if ((pa eq 0.) and (ra_c eq 0.) and (dec_c eq 0.)) then begin
	      print, 'Error: Wrong galaxy name, or No records for the galaxy !'
	      return
	  endif

	; get the position in [pixel] of galaxy center  
	  cube=readfits(infile,hdr)
	  ; noise=mad(cube,/finite)
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
	
	  ad2xy, ra_c, dec_c, astr, xc, yc
	  print,'===================================='
	  print,'Center of Galaxy in Pixel:', xc, yc
	  print,'===================================='

  endif


; READ ROTFILE - circular velocity at each radius
  readcol, rotdat, rad, vcirc, F='f, f', /silent  ; in [pc], [km/s]
  omega = vcirc/rad    ; km/s pc^-1
  deriv_omega = deriv(rad, omega)  
  deriv_omega2 = deriv(rad, omega^2.)  

  
; CALCULATE SHEAR PARAMETERS
  A_array = make_array(nclumps)
  T_array = make_array(nclumps)
  var_diff_array = make_array(nclumps)
  dvar_diff_array = make_array(nclumps)
  omega_array = make_array(nclumps)


  for i = 0, nclumps -1 do begin
    Rc = R[i]            ; cloud size in [pc]
    d0 = distance[i]     ; cloud distance in [pc]
    omega0 = INTERPOL(omega, rad, d0)
    A0 = -0.5 * d0 * INTERPOL(deriv_omega, rad, d0)
    B0 = A0 - omega0
    T0 = 4 * omega0 * A0
    ;T0 = -1. * d0 * INTERPOL(deriv_omega2, rad, d0)
	; Mlum0 = Mlum[i]
    ; G = 4.302e-3           ; in [pc Msun^-1 *(km/s)^2]


    this_mom = props[i]
    mom1x = this_mom.mom1x
    mom1y = this_mom.mom1y
    maj=(mom1x-xc)*sin(pa*!dtor)-(mom1y-yc)*cos(pa*!dtor) ; Eric: add (-1) here?
    min=(mom1x-xc)*cos(pa*!dtor)+(mom1y-yc)*sin(pa*!dtor) ; Eric: add (-1) here?
  	min=min/cos(inc*!dtor)
    this_c = maj/sqrt(maj^2.+min^2.)
    this_s = min/sqrt(maj^2.+min^2.)

    var_diff = (omega0^2*1/5.*Rc^2.*this_s^2+ $
		(omega0-2*A0)^2.*1/5.*Rc^2.*this_c^2) * sin(inc*!dtor)^2.

    dvar_diff = sqrt( (2*d_R[i]/Rc)^2 + (2*0.08)^2 ) * var_diff 
    ; inclination from Pandora's paper 13.6 +0.9 -1.2 deg. Fractional error of sin(i) is 8%

    
	; assign values to array
    A_array[i] = A0
	omega_array[i] = omega0
    T_array[i] = T0
	var_diff_array[i] = var_diff
	dvar_diff_array[i] = dvar_diff

  endfor

  openw, lun, outtable, /get_lun

  printf, lun, ';peaknum', ',', 'Omega', ',', 'T', ',','d(sig2)', ',', 'error', $
			format='(A10,a1,A5,a1,A1,a1,A7,a1,A7)'
  printf, lun, ' ', ',', 'km/s /pc', ',', '(km/s/pc)^2',',','km2/s2',',','km2/s2', $
			format='(A10,a1,A9,a1,A12,a1,A6,a1,A6)'
  for i = 0, nclumps-1 do begin
    printf, lun, ID[i],',', omega_array[i], ',', T_array[i], ',', var_diff_array[i], ',', dvar_diff_array[i], $
		format = '(I4,a1,F7.3,a1,F7.3,a1,F11.5,a1,F11.5)'
  endfor

  close, lun
  free_lun, lun


END

