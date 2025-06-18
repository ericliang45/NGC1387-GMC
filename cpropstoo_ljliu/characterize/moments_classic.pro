Function plane_function, x,y,p

; Z - velocity
; x - xaxis
; y - yaxis
;   p[0] - c
;   p[1] - a
;   p[2] - b

   ;logN=alog10((10.^x/p[0])^(p[1]+1.))
   Z = p[0] + p[1]*x + p[2]*y
   return, Z
End



function moments_classic $
   , props $
   , xin = x $
   , yin = y $
   , vin = v $
   , tin = t $
   , define_fields = define_fields $   
   , calculate = calculate $
   , moments = moments $
   , properties = properties $
   , hdr=hdr $
   , astr=astr $
   , vaxis=vaxis $
   , extrap = do_extrap $
   , extarg = targett $
   , fluxlin = fluxlin  $
   , do_bootstrap = do_bootstrap
  
;+
;
; Moment calculation module. Will either define fields, calculate
; moments, or convert moments into properties, depending on the
; calling syntax.
;
; This version implements the "classic" CPROPS prescription of
; Rosolowsky & Leroy (2006).
;
;-
;   20170213  LJ  add bootstrap uncertainties based on cprops_dist
;   20170213  LJ  add keyword do_bootstrap so do not calculate velocity
;				  gradient axis
;   20170214  LJ  add bootstrap uncertainties based on cprops_dist
;   20170223  LJ  add gaussian-fitted velocity dispersion which gets rid of
;				  large scale motions (more details see Utomo+ 2015)
;   20180108  LJ  add shape and convexity features 
;   20180307  LJ  smooth mom0 and mom1 by beamfwhm to fit rotation plane

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFINE FIELDS 
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; Add the necessary fields to the properties structure

; FOR CONVENIENCE, ASSIGN "not-a-number" TO THE nan VARIABLE
  nan = !values.f_nan

  if keyword_set(define_fields) then begin

;    ... ADD FIELDS FOR MOMENTS
     if keyword_set(moments) then begin
        
;       MOMENT 0  
        props = create_struct(props, "mom0", nan)
        props = create_struct(props, "mom0_extrap", nan)
        props = create_struct(props, "mom0_unit", "K*pix^3")
        
;       MOMENT 1
        props = create_struct(props, "mom1x", nan)
        props = create_struct(props, "mom1x_unit", "pix")
        
        props = create_struct(props, "mom1y", nan)
        props = create_struct(props, "mom1y_unit", "pix")

        props = create_struct(props, "mom1v", nan)
        props = create_struct(props, "mom1v_unit", "pix")
        
;       MOMENT 2
        props = create_struct(props, "mom2x", nan)
        props = create_struct(props, "mom2x_extrap", nan)
        props = create_struct(props, "mom2x_unit", "pix")
        
        props = create_struct(props, "mom2y", nan)
        props = create_struct(props, "mom2y_extrap", nan)
        props = create_struct(props, "mom2y_unit", "pix")
        
        props = create_struct(props, "mom2v", nan)
        props = create_struct(props, "mom2v_extrap", nan)
        props = create_struct(props, "mom2v_gauss", nan)
        props = create_struct(props, "mom2v_gauss_extrap", nan)
        props = create_struct(props, "mom2v_unit", "pix")

;       PRINCIPAL AXES CALCULATION
        props = create_struct(props, "mom2maj", nan)
        props = create_struct(props, "mom2maj_extrap", nan)
        props = create_struct(props, "mom2maj_unit", "pix")

        props = create_struct(props, "mom2min", nan)
        props = create_struct(props, "mom2min_extrap", nan)
        props = create_struct(props, "mom2min_unit", "pix")

        props = create_struct(props, "momposang", nan)
        props = create_struct(props, "momposang_unit", "rad")

        props = create_struct(props, "covar_xy", nan)
        
;       VELOCITY GRADIENT CALCULATION (ADDED DEC 4, 2014)
        props = create_struct(props, "velgrad", nan)
        props = create_struct(props, "velgrad_unit", "pix/pix")
        props = create_struct(props, "velposang", nan)
        props = create_struct(props, "velposang_unit", "rad")
        props = create_struct(props, "veldisp", nan)
        props = create_struct(props, "velchisq", nan)
        props = create_struct(props, "veldisp_unit", "pix")
        props = create_struct(props, "planefit_coef", [nan, nan, nan])
        props = create_struct(props, "planefit_coeferr", [nan, nan, nan])

;       SHAPE AND CONVEXITY FEATRURES
        props = create_struct(props, "texture", nan)
        props = create_struct(props, "Upar", nan)
        props = create_struct(props, "gradient", nan)
        props = create_struct(props, "surf_gradient", nan)
        props = create_struct(props, "shape", nan)
        props = create_struct(props, "convexity3d", nan)

;       DEFINE UNCERTAINTIES        
        props = create_struct(props, "radrms_extrap_deconv_uc", nan)
        props = create_struct(props, "vrms_extrap_deconv_uc", nan)  
        props = create_struct(props, "vrms_gauss_extrap_deconv_uc", nan)  
        props = create_struct(props, "lum_extrap_uc", nan) 
        props = create_struct(props, "mass_extrap_uc", nan)  
        props = create_struct(props, "virmass_extrap_deconv_uc", nan)
        props = create_struct(props, "virmass_gauss_extrap_deconv_uc", nan)

        return, props

     endif

;    ... ADD FIELDS FOR PROPERTIES
     if keyword_set(properties) then begin

;       FLUX-DERIVED QUANTITIES
        props = create_struct(props, "flux", nan)
        props = create_struct(props, "flux_extrap", nan)
        props = create_struct(props, "flux_unit", "K*km/s*as^2")

        props = create_struct(props, "lum", nan)
        props = create_struct(props, "lum_extrap", nan)
        props = create_struct(props, "lum_unit", "K*km/s*pc^2")

        props = create_struct(props, "mass", nan)
        props = create_struct(props, "mass_extrap", nan)
        props = create_struct(props, "mass_unit", "Msun")

;       POSITIONAL QUANTITIES
        props = create_struct(props, "xpos", nan)
        props = create_struct(props, "xpos_unit", "deg")

        props = create_struct(props, "ypos", nan)
        props = create_struct(props, "ypos_unit", "deg")

        props = create_struct(props, "vpos", nan)
        props = create_struct(props, "vpos_unit", "km/s")
  
;       SIZE QUANTITIES
        props = create_struct(props, "xrms", nan)
        props = create_struct(props, "xrms_extrap", nan)
        props = create_struct(props, "xrms_extrap_deconv", nan)
        props = create_struct(props, "xrms_unit", "pc")

        props = create_struct(props, "yrms", nan)
        props = create_struct(props, "yrms_extrap", nan)
        props = create_struct(props, "yrms_extrap_deconv", nan)
        props = create_struct(props, "yrms_unit", "pc")

        props = create_struct(props, "vrms", nan)
        props = create_struct(props, "vrms_deconv", nan)
        props = create_struct(props, "vrms_extrap", nan)
        props = create_struct(props, "vrms_extrap_deconv", nan)
        props = create_struct(props, "vrms_gauss", nan)
        props = create_struct(props, "vrms_gauss_deconv", nan)
        props = create_struct(props, "vrms_gauss_extrap", nan)
        props = create_struct(props, "vrms_gauss_extrap_deconv", nan)
        props = create_struct(props, "vrms_unit", "km/s")

;       PRINCIPAL AXIS QUANTITIES
        props = create_struct(props, "majrms", nan)
        props = create_struct(props, "majrms_extrap", nan)
        props = create_struct(props, "majrms_extrap_deconv", nan)
        props = create_struct(props, "majrms_unit", "pc")

        props = create_struct(props, "minrms", nan)
        props = create_struct(props, "minrms_extrap", nan)
        props = create_struct(props, "minrms_extrap_deconv", nan)
        props = create_struct(props, "minrms_unit", "pc")

        props = create_struct(props, "posang_extrap_deconv", nan)
        props = create_struct(props, "posang_unit", "rad")

        props = create_struct(props, "resolved_extrap", 0B)

;       PHYSICAL QUANTITIES
        props = create_struct(props, "radrms_extrap_deconv", nan)
        props = create_struct(props, "radrms_unit", "pc")

        props = create_struct(props, "virmass_extrap_deconv", nan)
        props = create_struct(props, "virmass_gauss_extrap_deconv", nan)
        props = create_struct(props, "virmass_unit", "Msun")

;       PHYSICAL UNCERTAINTIES

;       added by Eric

        props = create_struct(props, "gal_cen_x", nan)
        props = create_struct(props, "gal_cen_x_unit", "pixel")
        props = create_struct(props, "gal_cen_y", nan)
        props = create_struct(props, "gal_cen_y_unit", "pixel")

        props = create_struct(props, "gal_pa", nan)
        props = create_struct(props, "gal_pa_unit", "deg")
        props = create_struct(props, "gal_inc", nan)
        props = create_struct(props, "gal_inc_unit", "deg")

        props = create_struct(props, "gal_dv", nan)
        props = create_struct(props, "gal_dv_unit", "km/s")

        props = create_struct(props, "r_gal", nan)
        props = create_struct(props, "r_gal_unit", "pc")

        props = create_struct(props, "resolve_spatial", nan)
        props = create_struct(props, "resolve_spectral", nan)

        return, props

     endif

  endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; CALCULATE MOMENTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

;+ 
;
; Carry out the original CPROPS moment calculations.
;
;-

  if keyword_set(calculate) and keyword_set(moments) then begin

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; DEFAULTS & DEFINITIONS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

;    SET THE FUNCTIONAL FORM FOR THE FLUX EXTRAPOLATION
     if n_elements(fluxlin) eq 0 then $
        square = 1 $
     else $
        square = 0

;    SET THE FUNCTIONAL FORM FOR THE FLUX EXTRAPOLATION
     if keyword_set(do_extrap) then $
        do_extrap = 1 $
     else $
        do_extrap = 0

;    SET THE TARGET FOR THE CURVE OF GROWTH
     if n_elements(targett) eq 0 then $
        targett = 0
     
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; CUMULATIVE CALCULATIONS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; Calculate cumulative moments for the provided x, y, v, t vector.

;    ZEROTH MOMENT
     mom0t = total(t, cumul=do_extrap)

;    FIRST MOMENT
     mom1x = total(x*t, cumul=do_extrap)/total(t, cumul=do_extrap)
     mom1y = total(y*t, cumul=do_extrap)/total(t, cumul=do_extrap)
     mom1v = total(v*t, cumul=do_extrap)/total(t, cumul=do_extrap)

;    SECOND MOMENT
     term1x = total(t*double(x)^2., cumul=do_extrap)
     term2x = (total(t*double(x), cumul=do_extrap))^2./mom0t
     mom2x = sqrt((term1x - term2x)/mom0t)
     zeroind = where(abs(term1x - term2x) lt 1.d-10, num)
     if (num gt 0) then $
        mom2x[zeroind] = 0.0
     
     term1y = total(t*double(y)^2., cumul=do_extrap)
     term2y = (total(t*double(y), cumul=do_extrap))^2./mom0t
     mom2y = sqrt((term1y - term2y)/mom0t)
     zeroind = where(abs(term1y - term2y) lt 1.d-10, num)
     if (num gt 0) then $
        mom2y[zeroind] = 0.0  

     term1v = total(t*double(v)^2., cumul=do_extrap)
     term2v = (total(t*double(v), cumul=do_extrap))^2./mom0t
     mom2v = sqrt((term1v - term2v)/mom0t)
     zeroind = where(abs(term1v - term2v) lt 1.d-10, num)
     if (num gt 0) then $
        mom2v[zeroind] = 0.0

;    SAVE TO STRUCTURE
     props.mom0 = mom0t[n_elements(mom0t)-1]

     props.mom1x = mom1x[n_elements(mom1x)-1]
     props.mom1y = mom1y[n_elements(mom1y)-1]
     props.mom1v = mom1v[n_elements(mom1v)-1]

     props.mom2x = mom2x[n_elements(mom2x)-1]
     props.mom2y = mom2y[n_elements(mom2y)-1]
     props.mom2v = mom2v[n_elements(mom2v)-1]
     
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; PRINCIPAL AXES
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; Perform the principal axes decomposition and record those values.

     ellfit, x=x, y=y, wt=t, posang=posang

;    CALCULATE THE ROTATED X AND Y
     xrot = x*cos(posang)+y*sin(posang)
     yrot = -x*sin(posang)+y*cos(posang)    

;    FIRST MOMENT
     mom1xrot = total(xrot*t, cumul=do_extrap)/total(t, cumul=do_extrap)
     mom1yrot = total(yrot*t, cumul=do_extrap)/total(t, cumul=do_extrap)

;    SECOND MOMENT
     term1xrot = total(t*double(xrot)^2., cumul=do_extrap)
     term2xrot = (total(t*double(xrot), cumul=do_extrap))^2./mom0t
     mom2xrot = sqrt((term1xrot - term2xrot)/mom0t)
     zeroind = where(abs(term1xrot - term2xrot) lt 1.d-10, num)
     if (num gt 0) then $
        mom2xrot[zeroind] = 0.0
     
     term1yrot = total(t*double(yrot)^2., cumul=do_extrap)
     term2yrot = (total(t*double(yrot), cumul=do_extrap))^2./mom0t
     mom2yrot = sqrt((term1yrot - term2yrot)/mom0t)
     zeroind = where(abs(term1yrot - term2yrot) lt 1.d-10, num)
     if (num gt 0) then $
        mom2yrot[zeroind] = 0.0    
     
;    SAVE TO STRUCTURE
     props.momposang = posang     

     props.mom2maj = mom2xrot[n_elements(mom2xrot)-1]
     props.mom2min = mom2yrot[n_elements(mom2yrot)-1]

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; CURVE OF GROWTH EXTRAPOLATION
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; Perform the curve-of-growth extrapolation.

     if do_extrap then begin

        ex_mom2x = extrap(t, mom2x, targett = targett, /fast $
                          , scatter = e_mom2x, /weight)
        
        ex_mom2y = extrap(t, mom2y, targett = targett, /fast $
                          , scatter = e_mom2y, /weight)
        
        ex_mom2v = extrap(t, mom2v, targett = targett, /fast $
                          , scatter = e_mom2v, /weight)

        ex_mom2xrot = extrap(t, mom2xrot, targett = targett, /fast $
                             , scatter = e_mom2xrot, /weight)
        
        ex_mom2yrot = extrap(t, mom2yrot, targett = targett, /fast $
                             , scatter = e_mom2yrot, /weight)
        
        ex_mom0t = extrap(t, mom0t, targett = targett, /fast, square = square, $
                          scatter = e_mom0t, /weight)

        if (ex_mom0t lt max(mom0t, /NAN)) then begin
           help, calls = calls
           ; if n_elements(calls) eq 2 then begin
              ; print, 'ERROR in EXTRAPOLATION!!! Extrapolated flux is LESS THAN measured flux or infinite. That is BAD!'
              ; print, 'Defaulting to linear extrapolation...'
           ; endif
           ex_mom0t = extrap(t, mom0t, targett = targett, /fast)
        endif

;       ASSIGN TO STRUCTURE     
        props.mom0_extrap = ex_mom0t

        props.mom2x_extrap = ex_mom2x
        props.mom2y_extrap = ex_mom2y
        props.mom2v_extrap = ex_mom2v
		props.mom2v_gauss = (calc_rmsv_gauss(x, y, v, t))[1]
		props.mom2v_gauss_extrap = (calc_rmsv_gauss(x, y, v, t))[0]

        props.mom2maj_extrap = ex_mom2xrot
        props.mom2min_extrap = ex_mom2yrot

     endif

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; VELOCITY GRADIENT AXIS (added on Dec 4, 2014 by schruba@mpe.mpg.de)
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; Eric note:
; original code only ~60 lines, see https://github.com/akleroy/cpropstoo/blob/master/characterize/moments_classic.pro; 
; here, added to ~300 lines in total

   if keyword_set(do_bootstrap) then begin  ; added, don't run this routine for bootstrapping
     return, props ; added
   endif else begin ; added
	 pad = 2 ; added, increase "pad" parameter from 1 to 2

;      cubify $
;         , x=x, y=y, v=v, t=t $
;         , cube = minicube $
;         , pad = pad $ ; changed, see above
;         , location = location $
; 		, mask = mask $ ; added, to calculate geometry parameters such as convexity
;         , /silent

;      l3d = array_indices(minicube, location)

;      ; print, t[-1],x[-1],y[-1],v[-1]
;      ; print, minicube[l3d[0,-1],l3d[1,-1],l3d[2,-1]],l3d[0,-1],l3d[1,-1],l3d[2,-1]

;      rdhd, hdr, s = hdstruct ; added, to get vaxis in physical units, same for lines below
;      vaxis_old = hdstruct.v ; added
;      if strmatch(sxpar(hdr,'CUNIT3'),'*KM/S*') or $ ; added
;        strmatch(sxpar(hdr,'CUNIT3'),'*km/s*') then begin ; added
;        vaxis_old = vaxis_old*1e3    ; in [km/s] ; added, this is valid, since "rdhd" routine outputs velocity in mega-m/s, confusingly 
;      endif ; added

;      sz = size(minicube)
;      xaxis = (indgen(sz[1])+min(x)-pad) ; changed, "1" to variable "pad"
;      yaxis = (indgen(sz[2])+min(y)-pad) ; changed, "1" to variable "pad"
;      vaxis = (indgen(sz[3])+min(v)-pad) ; changed, "1" to variable "pad"
; 	 vaxis = vaxis*(vaxis_old[1]-vaxis_old[0])+vaxis_old[0]   ; changed, in [km/s]
;      xmap = (intarr(sz[2])+1) ## xaxis
;      ymap = yaxis ## (intarr(sz[1])+1)
;      vcube = rebin(reform(vaxis,1,1,sz[3]), sz[1:3])       ; in [km/s]

;      mom0 = total(minicube,3,/nan)
;      mom1 = total(vcube*minicube,3,/nan) / mom0     ; in [km/s]
; 	 mask2d = total(mask,3, /nan)/(total(mask, 3, /nan) > 1) ; added, to calculate geometry parameters such as convexity

; 	 covcube = minicube * !values.f_nan ; added, to produce "mom2_arr" (moment 2 map), similar to "sqrt(total(errcube,3,/nan))" in the original code; but both are wrong since (vcube-covcube) has values at many more voxels than the cloud actually has.
; 	 for i = 0, sz[3]-1 do begin  ; added, see above
; 	   covcube[*,*,i] = mom1        ; in [km/s] ; added, see above
; 	 endfor ; added, see above 
; 	 mom2_arr = sqrt(total(minicube*(vcube-covcube)^2, 3, /nan)/total(minicube, 3, /nan))     ; in [km/s] ; added, see above

	 
; 	 extast, hdr, astr ; added, to convolve mini-moment 1 map of clouds, obsolete now (only "beamfwhm_pix" still in use now)
; 	 xy2ad, [0,1], [0,0], astr, ra, dec  ; added, see above
; 	 degperpix = sphdist(ra[0], dec[0], ra[1], dec[1], /deg)  ; added, see above
; 	 beamfwhm_deg = sqrt(sxpar(hdr, "BMAJ") * sxpar(hdr, "BMIN"))  ; added, see above
; 	 beamfwhm_pix = beamfwhm_deg / degperpix  ; added, see above
;      psf=psf_gaussian(npixel=[sz[1],sz[2]], $  ; added, see above
;           FWHM = beamfwhm_pix, /normal)  ; added, see above
;      ;mom1=convol(mom1, psf, /nan, /edge_mirror, /center, $  ; added, see above
;      ;           missing=NaN, /normalize)  ; added, see above
;      ;mom0=convol(mom0, psf, /nan, /edge_mirror, /center, $  ; added, see above
;      ;           missing=NaN, /normalize)  ; added, see above

; ;    ESTIMATE UNCERTAINTY IN VELOCITY FIELD
;      errcube = (vcube - rebin(mom1, sz[1:3]))^2
;      mom1err = props.noise/mom0 * sqrt(total(errcube,3,/nan))


; ;    FIT A PLANE TO THE MOMENT MAP (REQ FREUDENREICH ROUTINES)
; ;    resid : cloud velocity dispersion after gradient correction
;      if (keyword_set(robust)) then begin   ; weird, no keyword called robust (comment from Eric)
;         fitmap = ROB_MAPFIT(mom1, 1, coef, resid)
;      endif else begin
;         mom1vec = reform(mom1[l3d[0,*],l3d[1,*]])
;         mom2vec = reform(mom2_arr[l3d[0,*],l3d[1,*]])  ; added
;         weights = 1.0/(mom1err^2)
;         wt = reform(weights[l3d[0,*],l3d[1,*]])

;         ind = WHERE_XYZ(finite(mom1), XIND=xind ,YIND=yind)  ; added by Eric, to get a vector for fitting (each LoS only sampled once)
;         mom1vec_new = mom1[ind]   ; added, see above
;         ; print, xind[0], yind[0], mom1vec_new[0], mom1[xind[0], yind[0]]  ; added, see above, for testing
;         ; print, xind[1], yind[1], mom1vec_new[1], mom1[xind[1], yind[1]]  ; added, see above, for testing
;         ; print, xind[2], yind[2], mom1vec_new[2], mom1[xind[2], yind[2]]  ; added, see above, for testing

;         ; if total(finite(wt)) ge 10 then begin ; original line ; removed/changed by Eric
;         if total(finite(mom1vec_new)) ge 10 then begin ; added by Eric

;            ;coef = PLANEFIT(x, y, mom1vec, wt, vfit)   ; removed/changed to mpfit2dfun
; 		   ;if n_elements(vfit) gt 0 then begin  ;  removed/changed to mpfit2dfun
;            ;  resid = stddev(vfit - mom1vec)  ; removed/changed to mpfit2dfun
; 		   ;endif else begin  ; removed/changed to mpfit2dfun
;            ;  coef = fltarr(3) * !values.f_nan  ; removed/changed to mpfit2dfun
;            ;  resid = !values.f_nan  ; removed/changed to mpfit2dfun
; 		   ;endelse

; 		   p0 = [1300,1.0,0] ; need to be updated for a more generic case
;            coef = mpfit2dfun('plane_function',x,y,mom1vec,mom2vec,p0,weights=wt,perror=perror,quiet=1) ; original
;            ; coef = mpfit2dfun('plane_function',x,y,mom1vec,1.0,p0,perror=perror,quiet=1) ; Eric
;            ; coef = mpfit2dfun('plane_function',xind,yind,mom1vec_new, FLTARR((size(xind))[1])+1.0 ,p0,perror=perror,quiet=1) ; Eric v2, error=1.0 m/s, half channel width
;            resid = !values.f_nan
;         endif else begin
;            coef = fltarr(3) * !values.f_nan
;            resid = !values.f_nan
;            perror = !values.f_nan
;         endelse
;      endelse

; ;    CALCULATE GRADIENT AND POSITION ANGLE
;      velgrad = sqrt(coef[1]^2+coef[2]^2) ; gradient in 'km/s/pixel' units
;      velposang = atan(coef[2],coef[1])   ; angle from x-axis (increasing R.A.)
;                                          ; to y-axis (increasing DEC)

; ;    All following are added by Lijie to calculated reduced chi2 of the fitting (in 1D dimension, v vs distance)

; ;    CALCULATE RESIDUATE CHISQ 
;      x0 = total(t*x, /NAN)/total(t, /NAN)
; 	 y0 = total(t*y, /NAN)/total(t, /NAN)

; 	 ; distance to rotation axis; d_v means distance to the kinematic axis
; 	 a = coef[1]/coef[2]
; 	 b = 1.
; 	 c = -(coef[1]/coef[2] * x0 + y0)
;   	 sign = (x-x0)*coef[1]+(y-y0)*coef[2]
;      sign = sign/abs(sign)
;      d_v = sign * abs(a*x + b*y + c)/sqrt(a^2. + b^2.)  ;in [pixel]

;      ; bin points
;      nbins = 10.
;      yhist=histogram(d_v, nbins = nbins, locations=xbins_hist, reverse_indice = R)
;      binsize=(max(d_v)-min(d_v))/(nbins-1)

; 	 xbins=[]
; 	 ybins=[]
; 	 xerr=[]
; 	 yerr=[]
; 	 prev=[]
; 	 xprev=[]
; 	 yprev=[]
; 	 the_total = n_elements(uniq(mom1vec[sort(mom1vec)]))
; 	 for i =0, nbins-1 do begin
; 	   if R[i] ne (R[i+1]) then begin
; 		 the_array = [mom1vec[R[R[i]:R[i+1]-1]],yprev]
; 	     if (stddev([mom1vec[R[R[i]:R[i+1]-1]],yprev]) gt 0.0) and $
; 			 (n_elements(uniq(the_array[sort(the_array)])) ge $
; 			 (2. > the_total/nbins)) then begin
; 	       xbins=[xbins, mean([d_v[R[R[i]:R[i+1]-1]], xprev])]
; 	       ybins=[ybins, mean([mom1vec[R[R[i]:R[i+1]-1]], yprev])]
; 	       xerr=[xerr, 0.0]
; 	       yerr=[yerr, stddev([mom1vec[R[R[i]:R[i+1]-1]], yprev])]
	 
; 	       xprev=[]
; 	       yprev=[]
; 	     endif else begin
; 		   xprev = [xprev, d_v[R[R[i]:R[i+1]-1]]]
; 	       yprev = [yprev, mom1vec[R[R[i]:R[i+1]-1]]]
; 	     endelse
; 	   endif
; 	 endfor

;      ; commented out by Eric
;      ; fit bin points and calculate reduced chisq
; 	 ; if n_elements(d_v) ge 2 then begin
;   ;       b_temp = mom1vec[UNIQ(mom1vec, SORT(mom1vec))] ; added by Eric
;   ;       if size(b_temp, /N_DIMENSIONS) gt 0 then begin
; 		;  sixlin, d_v, mom1vec, a1, siga1, b1, sigb1
;   ;       endif else begin
;   ;        a1 = !values.f_nan
;   ;        b1 = !values.f_nan
;   ;       endelse            
; 	 ; endif else begin
; 		;  a1 = !values.f_nan
; 		;  b1 = !values.f_nan
; 	 ; endelse
;          a1 = !values.f_nan
;          b1 = !values.f_nan

;         velchisq = !values.f_nan
;         ; commented by Eric
; 	 ; if n_elements(ybins) gt 2. then begin
; 		;  velchisq = total((ybins-(a1[1]+xbins*b1[1]))^2./yerr^2.) $
; 		; 	 /(n_elements(ybins)-2.)
; 	 ; endif else begin
; 		;  velchisq = !values.f_nan
; 	 ; endelse

; ;    SAVE RESULTS
;      props.velgrad = velgrad
;      props.velgrad_unit = 'kms/pix'   ; changed, since v has always been in physical units in this version
;      props.velposang = velposang
;      props.velposang_unit = 'rad'
;      props.veldisp = resid
;      props.veldisp_unit = 'kms'  ; changed, since v has always been in physical units in this version

;      props.velchisq = velchisq      ; added
; 	 props.planefit_coef = coef         ; kms/pix; added
; 	 props.planefit_coeferr = perror    ; kms/pix; added



; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; SHAPE AND CONVEXITY FEATURES (added on Jan 8, 2018 by ljliu.astro@gmail.com) 
; More details see Lin et al. (2003) - "A hybrid 3D watershed algorithm 
; incorporating gradient cues and object models for automatic segmentation of
; nuclei in confocal image stacks"
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    CALCULATE TEXTURE - the standard deviation of intensities of all pixels
;    inside the object
     t0 = mean(t,/nan)
	 texture = sqrt(total((t - t0)^2.)/(n_elements(t)-1.))

;;    CALCULATE THE SHAPE
     thepad = 2
     cubify $
     , x=x, y=y, v=v, t=t $
     , cube = theminicube $
     , pad = thepad $
     , location = thelocation $
     , mask = mask


     thesz = size(theminicube)
	 thexaxis = (indgen(thesz[1])+min(x)-thepad)
	 theyaxis = (indgen(thesz[2])+min(y)-thepad)
	 thexmap = (intarr(thesz[2])+1) ## thexaxis
	 theymap = theyaxis ## (intarr(thesz[1])+1)

;    CALCULATE THE Upar
	 mask1 = shift(mask, 1,0,0)
     mask2 = shift(mask, -1,0,0)
     mask3 = shift(mask, 0,1,0)
     mask4 = shift(mask, 0,-1,0)
     mask5 = shift(mask, 0,0,1)
     mask6 = shift(mask, 0,0,-1)
     mask_sum = mask1+ mask2 + mask3 + mask4 + $
     mask5 + mask6 + mask

	 inside = where(mask_sum eq 7, inside_ct)
	 Vpar = n_elements(thelocation)
	 Qpar = Vpar - inside_ct ; number of the boundary pixels
	 Upar = Qpar^3./(64.*!pi*Vpar^2.)

;    CALCULATE THE GRADIENT
     l3d = array_indices(theminicube, thelocation)
     grad_array = thelocation*0.0
     nan_ind= where(theminicube ne theminicube)
     theminicube1 = theminicube
     theminicube1[nan_ind] = 0.0
     for i = 0, n_elements(thelocation)-1 do begin
        x3d = l3d[0,i]
        y3d = l3d[1,i]
        z3d = l3d[2,i]
        term1 = (theminicube1[x3d+1,y3d,z3d]-theminicube1[x3d-1,y3d,z3d])/2.
        term2 = (theminicube1[x3d,y3d+1,z3d]-theminicube1[x3d,y3d-1,z3d])/2.
        term3 = (theminicube1[x3d,y3d,z3d+1]-theminicube1[x3d,y3d,z3d-1])/2.
        grad_array[i]= (term1^2. + term2^2. + term3^2.)/3.
     endfor
     gradient = sqrt(total(grad_array*t, /nan)/(total(t,/nan)))

;    CALCULATE THE SURFACE GRADIENT
     surf_ind = where((mask eq 1) and (mask_sum lt 7), surface_ct)
     surf_t = theminicube[surf_ind]
	 surf_t0 = mean(surf_t, /nan)
	 surf_gradient = sqrt(total((surf_t - surf_t0)^2.*surf_t, /nan)/ $
			(total(surf_t, /nan)))


;    CALCULATE 2D CONVEXITY
     mask2d = total(mask,3, /nan)/(total(mask, 3, /nan) > 1)
     locat2d = where(mask2d gt 0, area)
	 if area gt 0 then begin
	     x2d = thexmap[locat2d]
	     y2d = theymap[locat2d]
	 
	     convex_area = convexhull_area(x2d,y2d)
	     shape = area/convex_area
	 endif


;    CALCULATE 3D CONVEXITY
     themom0 = total(theminicube,3,/nan)
	 I2d = themom0[locat2d]
	 x2d = thexmap[locat2d]
	 y2d = theymap[locat2d]
	 ;I2d = re_order(I2d)
	 ;iVolume = total(I2d)  ; in [intensity_unit*pixel*pixel]
	 ;convex_volume = convexhull_volume(x2d, y2d, I2d)
	 ;convexity3d = iVolume/convex_volume
	 convexity3d = calc_convexity(x2d,y2d,I2d)

     props.texture = texture
	 props.Upar = Upar
	 props.gradient = gradient
	 props.surf_gradient = surf_gradient
     props.shape = shape
     props.convexity3d = convexity3d

	 ;print,'shape,convexity3d,texture, Upar',shape,convexity3d,texture, Upar



; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; BE DONE
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     return, props
   endelse   ; keyword do_bootstrap

  endif
  
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; CONVERT MOMENTS TO PROPERTIES
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

;+ 
;
; Convert from moments to properties. At this point we can assume a
; user-defined alpha and distance and header information, so that we
; know the astromery and velocity.
;
;-

  if keyword_set(calculate) and keyword_set(properties) then begin

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; FLUX-LIKE QUANTITIES
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     mom_to_flux = props.chanwidth_kms*(props.degperpix*3600.)^2
     props.flux = props.mom0*mom_to_flux
     props.flux_extrap = props.mom0_extrap*mom_to_flux

     mom_to_lum = props.chanwidth_kms*(props.pcperpix)^2
     props.lum = props.mom0*mom_to_lum
     props.lum_extrap = props.mom0_extrap*mom_to_lum

     mom_to_mass = props.chanwidth_kms*(props.pcperpix)^2*props.alpha
     props.mass = props.mom0*mom_to_mass
     props.mass_extrap = props.mom0_extrap*mom_to_mass

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; POSITIONAL QUANTITIES
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     xy2ad, props.mom1x, props.mom1y, astr, ra_mom1, dec_mom1
     props.xpos = ra_mom1
     props.ypos = dec_mom1

     channum = findgen(n_elements(vaxis))
     props.vpos = interpol(vaxis, channum, props.mom1v)

     ; added by Eric
    props.r_gal=calc_distance(x1=props.gal_cen_x, y1=props.gal_cen_y, x2=props.mom1x, y2=props.mom1y , $
            in_pa=props.gal_pa, in_inc=props.gal_inc, pcperpix=props.pcperpix)
  
     
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; SIZE-LIKE QUANTITIES
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     props.xrms = props.mom2x*props.pcperpix
     props.xrms_extrap = props.mom2x_extrap*props.pcperpix
     props.xrms_extrap_deconv = $
        sqrt(props.mom2x_extrap^2 - (props.beamfwhm_pix/props.sig_to_fwhm)^2) * $
        props.pcperpix
     
     props.yrms = props.mom2y*props.pcperpix
     props.yrms_extrap = props.mom2y_extrap*props.pcperpix
     props.yrms_extrap_deconv = $
        sqrt(props.mom2y_extrap^2 - (props.beamfwhm_pix/props.sig_to_fwhm)^2) * $
        props.pcperpix

     props.vrms = props.mom2v*props.chanwidth_kms
     props.vrms_deconv = $
        sqrt(props.mom2v^2 - (1.0*props.chantosig)^2)*props.chanwidth_kms
     props.vrms_extrap = props.mom2v_extrap*props.chanwidth_kms
     props.vrms_extrap_deconv = $
        sqrt(props.mom2v_extrap^2 - (1.0*props.chantosig)^2)*props.chanwidth_kms

	 props.vrms_gauss_extrap = props.mom2v_gauss_extrap*props.chanwidth_kms
	 props.vrms_gauss_extrap_deconv = $
        sqrt(props.mom2v_gauss_extrap^2 - (1.0*props.chantosig)^2)*props.chanwidth_kms
	 props.vrms_gauss = props.mom2v_gauss*props.chanwidth_kms
	 props.vrms_gauss_deconv = $
        sqrt(props.mom2v_gauss^2 - (1.0*props.chantosig)^2)*props.chanwidth_kms

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; PRINCIPAL AXIS QUANTITIES
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     props.majrms = props.mom2maj*props.pcperpix
     props.majrms_extrap = props.mom2maj_extrap*props.pcperpix

     props.minrms = props.mom2min*props.pcperpix
     props.minrms_extrap = props.mom2min_extrap*props.pcperpix

;    DO THE TWO-D DECONVOLUTION (FACTORS MATCH SIGMA TO FWHM)
     deconvolve_gauss $
        , meas_maj = props.majrms_extrap*props.sig_to_fwhm $
        , meas_min = props.minrms_extrap*props.sig_to_fwhm $
        , meas_pa = props.momposang/!dtor $ ; INPUT AS DEGREES
        , beam_maj = props.bmaj_deg*!dtor*props.dist_pc $
        , beam_min = props.bmin_deg*!dtor*props.dist_pc $
        , beam_pa = props.bpa_deg $
        , src_maj = src_maj $
        , src_min = src_min $
        , src_pa = src_pa $
        , worked = worked $
        , point = point
     src_pa *= !dtor            ; SAVE AS RADIANS

     if worked then begin
        props.majrms_extrap_deconv = src_maj/props.sig_to_fwhm
        props.minrms_extrap_deconv = src_min/props.sig_to_fwhm
        props.posang_extrap_deconv = src_pa
        props.resolved_extrap = 1B
     endif else begin
        props.resolved_extrap = 0B
     endelse

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; VELOCITY GRADIENT AXIS IN PHYSICAL UNITS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   ; ;props.velgrad *= props.chanwidth_kms / props.pcperpix
   ; props.velgrad *= 1. / props.pcperpix  ; changed, since v has always been in physical units in this version
   ; props.planefit_coef *= 1./props.pcperpix        ; kms/pc ; added
   ; props.planefit_coeferr *= 1./props.pcperpix     ; kms/pc ; added
   ; ;props.veldisp *= props.chanwidth_kms  ; removed, since v has always been in physical units in this version
   ; props.velgrad_unit = 'km/s/pc'
   ; props.veldisp_unit = 'km/s'

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; PHYSICAL QUANTITIES
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     props.radrms_extrap_deconv = $
        props.rmstorad*sqrt(props.majrms_extrap_deconv* $
                            props.minrms_extrap_deconv)

     props.virmass_extrap_deconv = $
        props.vircoeff*props.radrms_extrap_deconv* $
        props.vrms_extrap_deconv^2 ;props.vrms_gauss_extrap_deconv^2 ; corrected by Eric

    ; added by Eric
    props.resolve_spatial = float(props.radrms_extrap_deconv ge props.beamfwhm_pc/2.)
    props.resolve_spectral = float(props.vrms_extrap_deconv ge props.gal_dv/2.)
     
     return, props

  endif

; SHOULD NOT GET HERE
  message, "No calculation carried out. Check calling sequence.", /info
  return, nan

end
   
