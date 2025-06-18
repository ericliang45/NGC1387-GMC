function calc_jtok, hdr=hdr $
                    , infile = infile $
                    , bmaj=bmaj $
                    , bmin=bmin $
                    , restfreq=restfreq $
                    , aips=aips $
                    , pixels_per_beam = pixels_per_beam


;+
;
; NAME:
;
;   CALC_JTOK()
;
; PURPOSE:
;   Helper function to turn cube in units of Jy/Beam into units of K.      
;
;
; CALLING SEQUENCE:
;   jtok = calc_jtok(hdr=hdr [,bmaj=bmaj, bmin=bmin,restfreq=restfreq,aips=aips]) 
;
; INPUTS:
;    HDR -- .fits header as an array of strings.
;    BMAJ -- (optional) Specify beam major axis (required without
;            header).
;    BMIN -- (optional) Specify beam minor axis (required without
;            header). 
;    restfreq -- (optional) The restfrequency of the observation
;                (required without header).  
;             
; KEYWORD PARAMETERS:
;    AIPS -- USE AIPS HEADER PARAMETERS
;
; OUTPUTS:
;   JTOK -- The conversion from Jansky/Beam to Kelvin
; MODIFICATION HISTORY:
;
;       Documentated -- Mon Nov 25  2013  Stephen Pardy 
;                     <spardy@astro.wisc.edu>
; 
;-


  ; PHYSICAL CONSTANTS
  c = 2.99792458d10             ; speed of light CGS
  h = 6.6260755d-27             ; Planck's constant CGS
  kb = 1.380658d-16             ; Boltzmann's const CGS

  if n_elements(hdr) gt 0 then begin
     if keyword_set(aips) then begin
        getaipsbeam, hdr, bmaj=bmaj, bmin=bmin
     endif else begin
	    if (sxpar(hdr, 'BMAJ') eq 0) or (sxpar(hdr, 'BMIN') eq 0) then begin
          if sxpar(hdr,'CASAMBM') then begin
			  print,'infile',infile
              beamtab=mrdfits(infile,1)
              if n_elements(beamtab) gt 1 then begin
				sxaddpar,hdr,'BMAJ', mean(beamtab.bmaj)/3600.
				sxaddpar,hdr,'BMIN', mean(beamtab.bmin)/3600.
				sxaddpar,hdr,'BPA', mean(beamtab.bpa)
			  endif
		  endif
		endif     ; end sxpar(hdr, 'BMAJ')
        if n_elements(bmaj) eq 0 then $
           bmaj = sxpar(hdr, 'BMAJ')
        if n_elements(bmin) eq 0 then $
           bmin = sxpar(hdr, 'BMIN')
     endelse
     if n_elements(restfreq) eq 0 then begin
        restfreq = sxpar(hdr, 'RESTFRQ', count=found_rf)
        if found_rf eq 0 then $
           restfreq = sxpar(hdr, 'RESTFREQ', count=found_rf)
        if found_rf eq 0 then $
           message, 'NO REST FREQUENCY!', /info
     endif
     getrot,hdr,rotation,cdelt
     pixels_per_beam = 2*!pi*bmaj*bmin/(8*alog(2)) / abs(cdelt[0]*cdelt[1])
  endif

  res_deg = sqrt(bmaj*bmin)
  beaminster = !pi*(res_deg*!dtor/2.)^2/alog(2.)
  jtok = c^2/beaminster/1.d23/(2.*kb*restfreq^2.)

  ; 'jypb2k' from cpropsdist
  ;lam_mm = 3d11/restfreq
  ;Conversion Factors for the map to change from Jy/beam to K.
  ;jypb2k = 14*(lam_mm/(sqrt(bmaj*3600.*bmin*3600.)))^2
  ;print,'jypb2k',jypb2k

  ;print, '================'
  ;print, 'jtok', jtok
  ;print, '================'

  return, jtok

end
