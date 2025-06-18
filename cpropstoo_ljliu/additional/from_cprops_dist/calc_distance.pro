Function calc_distance, x1=x1, y1=y1, x2=x2, y2=y2, $
                  in_pa=pa, in_inc=inc, pcperpix=pcperpix
                  

;+
; NAME: 
;   calc_distance
;
; PURPOSE:
;   calculate the physical distance between two pixels in a disk galaxy
;   
;
; INPUTS:
;   x1, y1     -- the position in [pixel] of the first point (galaxy center)
;   x2, y2     -- the position in [pixel] of the second point
;   in_pa         -- the position angle in [degree] of major axis of galaxy disk 
;				  from the North to the East: PA = 0 corresponds to a galaxy
;			      whose longest axis is oriented North-South.
;                 But note that the east direction always points towards
;                 decreasing xPos direction
;   in_inc         -- in_inclination of disk in [degree]
;   pcperpix   -- pc per pixel
;    
; OUTPUTS:
;   the physical distance between [x1,y1] and [x2,y2] in galaxy disk
;
; 
; HISTORY:
;   20161214  LJ   introduced


; 1. maj, min
  in_pa = pa*!dtor
  in_inc = inc*!dtor


  maj=(x2-x1)*sin(in_pa)-(y2-y1)*cos(in_pa)
  min=(x2-x1)*cos(in_pa)+(y2-y1)*sin(in_pa)
  axis_ratio = cos(in_inc)
  min=min/axis_ratio

  if n_elements(pcperpix) ge 1 then begin
	  distance=sqrt((maj^2.)+(min^2.))*pcperpix   ; in [pc]
  endif else begin
	  distance=sqrt((maj^2.)+(min^2.))			  ; in [pixel]
  endelse

  return, distance

End
