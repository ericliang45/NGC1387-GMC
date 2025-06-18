Function re_order, in_array

; Generate an array whose value equals to the reverse order of its unique value
; Input:
;   array  - 1D array

  array_uniq = in_array[uniq(in_array, sort(in_array))]
  max = n_elements(array_uniq)
  
  out_array = in_array * !values.f_nan
  for i = 0, max-1 do begin
	  ind = where(in_array eq array_uniq[i])
	  cal = where(array_uniq gt array_uniq[i], ct)
	  out_array[ind] =  max - ct
  endfor


  return, out_array
  



END
