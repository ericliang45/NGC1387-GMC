Function gaussian_function, x, p

;  p - parameters of gaussian function
;    p[0] - normalization value
;    p[1] - mean value 
;    p[2] - rms dispersion

   z = (x-p[1])/p[2]
   yfit = p[0]*exp(-z^2./2.)

   return, yfit

END

function calc_rmsv_gauss, xin, yin, vin, tin


;+
; NAME:
;   CALC_RMSV_GAUSS
;
; PURPOSE:
;   To measure the cloud's velocity dispersions by Gaussian fitting
;   to the composite spectrum of each cloud through the following steps.
;	First, we calculate the offset of the mean velocity at all lines
;	of sight with the cloud, with respect to the mean velocity at the
; 	center of the cloud. Then, we shift each line of sight velocity spectrum 
;	to match the mean velocity of the central position of the cloud. 
;	Finally, we fit the composite spectrum (the average velocity profile
;   from each line of sight) with a Gaussian, the derived standard deviation
;	of the Gaussian fit is therefore the velocity despersion we need.
;	More details see Section 3.1 of Utomo et al. 2015 paper.
;
; CALLING SEQUENCE:
;   rmsv_gauss = calc_rmsv_gauss(xin, yin, vin, tin)
;
; INPUTS:
;	xin, yin, vin, tin -- Vectors containing the data points of cloud
;
; OUTPUT:
;   rmsv_gauss -- the gaussian-fitted velocity dispersion in pixel
;
; HISTORY:
;   20170221     LJ   introduced

;  Gengerate a cube for the cloud: spec_cube[x=0,1,...,y=0,1...,v=0,1,...]=t, 
;  vel_cube[x,y,v] = vel in [pixel] of original vaxis
   pad = 1
   cubify, x = xin, y = yin, v = vin, t = tin, cube=spec_cube, $
       mask = mask, pad = pad, location=location, /silent
   cubify, x = xin, y = yin, v = vin, t = vin, cube=vel_cube, $
       pad = pad, /silent
   mask2d=total(mask,3,/nan)

;  Center pixel of the cloud
   ;xctr0, yctr0 - the center pixel in original data cube
   xctr0 = floor(total(tin*xin, /NAN)/total(tin, /NAN))
   yctr0 = floor(total(tin*yin, /NAN)/total(tin, /NAN))
   ;xctr, yctr - the center pixel in new data cube (spec_cube, vel_cube)
   xctr = xctr0 - (min(xin) - pad)
   yctr = yctr0 - (min(yin) - pad)

;  Mean velocity at all lines of sight with cube
   sz=(size(spec_cube))[1:3]
   v0=make_array(sz[0],sz[1])*!values.f_nan

   for i = 0, sz[0]-1 do begin
	   for j = 0, sz[1]-1 do begin
         if mask2d[i,j] gt 0 then begin
		   ti = spec_cube[i,j,*]
		   vi = vel_cube[i,j,*]
		   v0[i,j] = floor(total(vi*ti, /nan)/total(ti, /nan))
         endif
	   endfor
   endfor

   mom0 = total(spec_cube, 3, /nan)
   v0 = total(vel_cube*spec_cube, 3, /nan) /mom0
   
   ; mean velocity at cloud center
   v0ctr = v0[xctr, yctr]
   if finite(v0ctr) eq 0 then $
		v0ctr=total(tin*vin, /NAN)/total(tin, /NAN)

;  PLANE FIT
   l3d = array_indices(spec_cube, location)
   mom1vec = reform(v0[l3d[0,*],l3d[1,*]])
   wt = xin*0.1
   ; coef = PLANEFIT(xin, yin, mom1vec, wt, vfit) ; commented by Eric


;  Shift each line of sight velocity spectrum to match the mean 
;  velocity of the central position of the cloud
   for i = 0, sz[0]-1 do begin
	   for j = 0, sz[1]-1 do begin
         if mask2d[i,j] gt 0 then begin
		   index = where(finite(vel_cube[i,j,*]))
           vel_cube[i,j,index] = vel_cube[i,j,index] - (v0[i,j]-v0ctr)
         endif
       endfor
   endfor

;  Do not do gaussian fit
   shift_v = vel_cube[location]
   mom0t = total(tin, cumul=1)
   term1v = total(tin*double(shift_v)^2., cumul=1)
   term2v = (total(tin*double(shift_v), cumul=1))^2./mom0t
   shift_mom2v = sqrt((term1v - term2v)/mom0t)
   zeroind = where(abs(term1v - term2v) lt 1.d-10, num)
   if (num gt 0) then $
      shift_mom2v[zeroind] = 0.0

   ex_shift_mom2v = extrap(tin, shift_mom2v, targett = 0, /fast $
                  , scatter = e_shift_mom2v, /weight)
   shift_mom2v_extrap = ex_shift_mom2v

   rms = shift_mom2v_extrap

   return, [rms, shift_mom2v[-1]]



;  Make a composite spectrum of velocity profile from each line of sight
   n_chan = max(vel_cube,/nan) - min(vel_cube,/nan) + 10
   channels = findgen(n_chan) + ((min(vel_cube,/nan) - 5) > 0)
   spectrum = make_array(n_chan)

   for i = 0, sz[0]-1 do begin
	   for j = 0, sz[1]-1 do begin
         if mask2d[i,j] gt 0 then begin
           this_spec = make_array(n_chan)
		   fin_index = where(finite(vel_cube[i,j,*]))
           thechans = value_locate(channels, vel_cube[i,j,fin_index])
           this_spec[thechans] = spec_cube[i,j,fin_index]
           spectrum = spectrum + this_spec
         endif
	   endfor
   endfor

   spectrum = spectrum / total((mask2d gt 0))

;  Gaussian fit to the composite spectrum
   ;yfit = gaussfit(channels, spectrum, A, nterms=3, sigma=sigma)
   ;rms = A[2]
   ;rms_uc = sigma[2]

   spec_index = where(abs(spectrum) gt 0, nchan_spec)

   expr = 'GAUSS1(x,p[0:2])'  ; p= [MEAN, SIGMA, AREA]
   yerr = make_array(n_elements(channels), value=1.)
   start = [v0ctr, 2, 2]
   pi = replicate({fixed:0, limited:[0,0], limits:[0.D, 0.D]},3)
   pi[1].limited=[1,1]
   chantosig = 1.0/sqrt(2.0*!pi)
   pi[1].limits=[chantosig, n_elements(channels)/2.]
   
   if (nchan_spec le 4) then begin
      pi[0].fixed = 1
   endif

   yfit = mpfitexpr(expr, channels, spectrum, yerr, start, $ 
	   parinfo = pi, bestnorm = bestnorm, perror = perror, $
	   status = status, errmsg = errmsg, /quiet)

   if status le 0 then begin
     print, 'Gaussian fit to the cloud spectra failed!'
     print, errmsg
     return, [rms, rms_uc]
   endif
 
   rms = yfit[1]
   rms_uc = perror[1]

 

   return, [rms, rms_uc]
   
   
end
   

