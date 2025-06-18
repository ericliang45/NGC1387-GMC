FUNCTION calc_convexity, x, y, t


; PURPOSE:
;   Calculate the convexity of a mom0 3D distribution


; INPUT:
;  x, y, z - input x, y, z array of a 3D-structure


; 1. convexity_1
  ivolume1 = total(t)
  convex_volume1 = convexhull_volume(x, y, t)
  convexity_1 = ivolume1/convex_volume1

; 2. convexity_2
  ; construct structure containing 80% flux
  ;conv_level = max(themom0, /nan) * 0.2
  ;conv_reg = label_region((themom0 ge conv_level),/all_neighbors,/ULONG)
  ;conv_max_ind = where(themom0 eq max(themom0))
  ;conv_ind = where(conv_reg[locat2d] eq conv_reg[conv_max_ind[0]])

  conv_ind = where(t ge max(t,/nan)*0.2)
  ivolume2 = total(t[conv_ind])
  convex_volume2 = convexhull_volume(x[conv_ind], y[conv_ind], t[conv_ind])
  convexity_2 = ivolume2/convex_volume2


  convexity = convexity_1 > convexity_2


  return, convexity

END
