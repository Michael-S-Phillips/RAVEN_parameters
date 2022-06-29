function mro_crism_lookupwv, wv, wvt, w=w

;function lookupwv, wv, det, lut, w=w


;-----------------------------------------------------------------------
; ABSTRACT
; function to find correct row for target wavelength
;
; INPUTS
; wv = target wavelength (nanometers)
; det = detector (0 = VNIR  1 = IR)
; lut = look up table array (3 by x) (Detector row, VNIR wvln, IR wvln)
; w = toggle return of wavelength rather than row.
;
; OUTPUT
; Returns row number closest to target wavelength
;
; HISTORY
; created 20006/04/07 by NRI for CRISM Summary Product support
;
; Note that wavelengths outside the min/max of the detector range will
; return row numbers for the closest wavelength in the detector.
;-----------------------------------------------------------------------


; Subtract wavelength from the values in the detector.
; Take absolute value, find minimum.
; That's the closest row to the target wavelength in the detector.

;case det of
;	0 : row = fix(lut(0,where ( abs(lut(1,*)-wv) eq min(abs(lut(1,*)-wv)))))
;	1 : row = fix(lut(0,where ( abs(lut(2,*)-wv) eq min(abs(lut(2,*)-wv)))))
;endcase


row = where( abs(wvt - wv) EQ min(abs(wvt - wv)) )

if (keyword_set(w)) then out = wvt[row] else out = row

return, out
end