function contour_prop, xin, yin, vin, tin, kernel, levels = levels_in $
  , noextrap = noextrap, all_neighbors = all_neighbors
;+
; NAME:
;   CONTOUR_PROP
;
; PURPOSE:
;   To calculate the moments of emission around a fixed local maximum.
;
; CALLING SEQUENCE:
;   moments = CONTOUR_PROP(X,Y,V,T, KERNELS [levels = levels, noextrap
;   = noextrap])
;
; INPUTS:
;   X,Y,V,T -- Vectors containing the data (from VECTORIFY.pro)
;   KERNELS -- Indices of the two kernels (in X, Y, V, T)
;
; KEYWORD PARAMETERS:
;   LEVELS -- The levels with which to contour the data.
;   NOEXTRAP -- Do not extrapolate moments to their values at 0K
;
; OUTPUTS:
;   MOMENTS -- An array of moments with the same number of elements as
;              KERNEL with tags definied by  CLOUDMOM.pro
;
; MODIFICATION HISTORY:
;
;	Fri Mar 24 11:36:22 2006, Erik 
;       Added all_neighbors keyword.
;
;       Documented -- Fri Sep 2 16:29:02 2005, Erik Rosolowsky
;                     <erosolow@asgard.cfa.harvard.edu>
;
;		
;
;-


  if n_elements(levels_in) eq 0 then levels = contour_values(tin) else $
    levels = levels_in

  nlmax = n_elements(lmaxin)
  x0 = xin[kernel] & y0 = yin[kernel] & v0 = vin[kernel]  

;   TO AVOID EXCESS WORK, IF WE HAVE A MINIMUM CONTOUR LEVEL, DROP
;   DATA BELOW THAT AT THE BEGINNING OF THE LOOP.

  contours = contourcloud(xin, yin, vin, tin, x0 = x0, y0 = y0, $
                          v0 = v0, clev = levels, all_neighbors = all_neighbors)
  use_contours = levels

;   MEASURE THE PROPERTIES AT EACH CONTOUR, THIS CAN BE QUITE TIME
;   INTENSIVE. THE OUTPUT IS AN ARRAY OF MOMENTS CALLED "momra." DO
;   THIS BY CALLING THE "cloudmom" ROUTINE TO CALCULATE THE CUMULATIVE
;   MOMENTS AND THEN EXTRAPOLATE THEM TO THE TARGET CONTOUR, "targett"
  mom = {rmsx:!values.d_nan, rmsy:!values.d_nan, $
         rmsv:!values.d_nan, flux:!values.d_nan, $
         ermsx:!values.d_nan, ermsy:!values.d_nan, $
         ermsv:!values.d_nan, eflux:!values.d_nan, $
         covar:!values.d_nan, number:0L, $
		 rmsv_gauss: !values.d_nan, $
		 rmsv_gauss_uc: !values.d_nan }   ;(LJ)

  momra = replicate(mom, n_elements(levels))

  for j = 0, n_elements(use_contours) - 1 do begin
    useind = where(contours ge use_contours[j], num)
    if num eq 0 then continue
    xuse = xin[useind]
    yuse = yin[useind]
    vuse = vin[useind]
    tuse = tin[useind]      
    
;   CALL THE "cloudmom" ROUTINE TO CALCULATE THE CUMULATIVE MOMENTS
;   AND THEN EXTRAPOLATE THEM TO THE TARGET CONTOUR, "targett"
    mom = cloudmom(xuse, yuse, vuse, tuse, mom0t = mom0t $
                   , targett = targett, noextrap = noextrap)
    momra[j] = mom
  endfor
  return, momra
end
