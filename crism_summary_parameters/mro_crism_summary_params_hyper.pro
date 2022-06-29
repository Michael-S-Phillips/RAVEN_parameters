;ABSTRACT
; Function to calculate summary products from CRISM hyperspectral data.
;
;INPUTS
;   cube = input cube -> dark subtracted, I/F cube or proxy.
;           X=along track (640),y = number of frames, z= #wavelengths
;           so, a given wvln slice of the cube is (*,*,wvln)
;   det  = detector flag (0=vnir, 1=ir, 2=joined)
;   wvt = wavelength array passed from ENVI header
;
; KEYWROORDS 
;	index_vec	Vector mapping to full list of possible parameters for the 
;	            specified detector, element =1 for desired params, =0 for
;	            params to be skipped.
;	ignore_val	Data ignore value; not used yet
;	band_names	A named variable for return of the band name list. If supplied,
;	            routine sets that list and returns without computing parameters.
;
;HISTORY
; May 06,2008  First version: MFM adapted from code supplied by B. Ehlmann 
; Jan 31,2009  FM: Added OLINDEX2 (improved olivine; Mark Salvatore)
; Jan 31,2009  FM: Added BD1980 (sulfite; Leah Roach)
; Jan 31,2009  FM: Added Doub2200 (2.22 & 2.28 um doublet; Leah Roach)
; Feb 09,2009  FM: Added BD2230 (2.23 micron band depth; Kim Lichtenberg)
; Feb 09,2009  FM: In BD2230 its R2198 not R21980; and init paramidx=NUM_VNIR_PARAM
;                  in IR section, instead of 0.
;
;-------------------------------------------------------------------------


function mro_crism_summary_params_hyper, cube, det, wvt, index_vec=index_vec, $
		ignore_val=ignore_val, band_names=band_names


; Update these constants if parameters are added
NUM_VNIR_PARAM = 0
NUM_IR_PARAM = 7
NUM_JOIN_PARAM = NUM_VNIR_PARAM + NUM_IR_PARAM


; There are no VNIR hyperspectral summary params defined yet. 
if (det eq 0) then begin
	; This should not ever happen; should be trapped in _event.pro
	print,' NO VNIR hyperspectral parameters defined.'
	return, 0
endif

if (~(keyword_set(index_vec))) then begin
	;if (det EQ 0) then index_vec = replicate(1, NUM_VNIR_PARAM)   ;S
    if (det EQ 1) then index_vec = replicate(1, NUM_IR_PARAM)     ; IRt
    if (det EQ 2) then index_vec = replicate(1, NUM_JOIN_PARAM)     ; joined
endif

if (keyword_set(band_names)) then begin
	; Return the band names list
    case(det) of
        0: begin
            band_names = ['']
        end
	
        1:begin
            band_names = ['OLINDEX2', 'BD1900R', 'BD1980', 'BD2200', 'Doub2200', 'BD2230', 'BD2500']
        end
	
        2:begin
            band_names = ['OLINDEX2', 'BD1900R', 'BD1980', 'BD2200', 'Doub2200', 'BD2230', 'BD2500']
        end    
    	else:begin
    	   band_names =  -1
        endelse
    endcase

    return_band_names = band_names[where(index_vec EQ 1)]
    return, return_band_names
endif


if (n_elements(ignore_val) eq 0) then ignore_val=65535.



sz = size(cube)  ; x = spatial (640 detector columns, unless binned)
		 		 ; y = spatial (# of frames along track)
		 		 ; z = spectral # wavelengths

NB = sz[3]       ; spectral bands
NX = sz[1]	 ; spatial detector
NY = sz[2]	 ; spatial along track


; For joined data:
if (det eq 2) then pcubej = make_array(nx, ny, NUM_JOIN_PARAM, value = 65535.)

; Make a copy of the input index vector:
index_vec_in = index_vec



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PARAMETERS USING VNIR WAVELENGTHS:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; None yet
;if ((det eq 0) or (det eq 2)) then begin
	; pcube = make_array(nx, ny, NUM_VNIR_PARAM, value = 65535.)

	; ---   future vnir parameters here   ----

	; if (det eq 2) then pcubej(*,*,0:10)=pcube

;endif                ; VNIR or J case



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PARAMETERS USING IR WAVELENGTHS:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (det eq 1) or (det eq 2) then begin

    ; Reset index_vec to IR params:
    if (det eq 2) then index_vec = index_vec_in[NUM_VNIR_PARAM:*]

	pcube = make_array(nx, ny, NUM_IR_PARAM, value = 65535.)


	;------------------------------------------------------------------------
	;Find the Improved Olivine Index OLINDEX2
	;------------------------------------------------------------------------
	paramidx = NUM_VNIR_PARAM
	;------------------------------------------------------------------------
	;Better identifies the 1-micron olivine absorption by accounting for
	;spectral slope.  Also corrects for the false identification of olivine
	;due to high albedo.  Incorrectly maps hydrated regions as olivine-rich
	;due to the associated decrease in reflectance beyond 2-microns.  Also,
	;like OLINDEX, this parameter falsely identifies pyroxenes as olivine.
	;future work will be done to correct both the false identifications of
	;hydration and pyroxene as olivine-rich.
	;------------------------------------------------------------------------
	if (index_vec[paramidx] EQ 1) then begin

		R2404 = mro_crism_lookupwv(2404,wvt)
		R2397 = mro_crism_lookupwv(2397,wvt)
		R2410 = mro_crism_lookupwv(2410,wvt)

		R1750 = mro_crism_lookupwv(1750,wvt)
		R1744 = mro_crism_lookupwv(1744,wvt)
		R1757 = mro_crism_lookupwv(1757,wvt)

		R1054 = mro_crism_lookupwv(1054,wvt)
		R1047 = mro_crism_lookupwv(1047,wvt)
		R1060 = mro_crism_lookupwv(1060,wvt)

		R1211 = mro_crism_lookupwv(1211,wvt)
		R1204 = mro_crism_lookupwv(1204,wvt)
		R1218 = mro_crism_lookupwv(1218,wvt)

		R1329 = mro_crism_lookupwv(1329,wvt)
		R1326 = mro_crism_lookupwv(1326,wvt)
		R1336 = mro_crism_lookupwv(1336,wvt)

		R1474 = mro_crism_lookupwv(1474,wvt)
		R1467 = mro_crism_lookupwv(1467,wvt)
		R1480 = mro_crism_lookupwv(1480,wvt)

		;Band averages
		AVG2404 = ((cube(*,*,R2404)+cube(*,*,R2397)+cube(*,*,R2410))/3)
		AVG1750 = ((cube(*,*,R1750)+cube(*,*,R1744)+cube(*,*,R1757))/3)
		AVG1054 = ((cube(*,*,R1054)+cube(*,*,R1047)+cube(*,*,R1060))/3)
		AVG1211 = ((cube(*,*,R1211)+cube(*,*,R1204)+cube(*,*,R1218))/3)
		AVG1329 = ((cube(*,*,R1329)+cube(*,*,R1326)+cube(*,*,R1336))/3)
 		AVG1474 = ((cube(*,*,R1474)+cube(*,*,R1467)+cube(*,*,R1480))/3)

		;Slope of the continuum
		contm = ((AVG2404-AVG1750)/(2.404-1.750))

		;y-intercept of the continuum
		yint = (AVG2404-(2.404*contm))

		;Expected values at 1.054, 1.211, 1.329, and 1.474 microns, respectively
		ex1054 = (contm*1.054)+yint
		ex1211 = (contm*1.211)+yint
		ex1329 = (contm*1.329)+yint
		ex1474 = (contm*1.474)+yint

		;Calculated band depths at 1.054, 1.211, 1.329, and 1.474 microns, respectively, with weights
		pcube(*,*,paramidx) = $
      (((ex1054-AVG1054)/ex1054)*0.1)+$
      (((ex1211-AVG1211)/ex1211)*0.1)+$
      (((ex1329-AVG1329)/ex1329)*0.4)+$
      (((ex1474-AVG1474)/ex1474)*0.4)

	endif


	;------------------------------------------------------------------------
	;Find the 1.90 micron H2O band depth (BD1900R): (wvs 1908-1941 (in), 1862-1875, 2112-2126)
	;------------------------------------------------------------------------
	paramidx++
	if (index_vec[paramidx] EQ 1) then begin
		R1908 = mro_crism_lookupwv(1908,wvt)
		R1914 = mro_crism_lookupwv(1914,wvt)
		R1921 = mro_crism_lookupwv(1921,wvt)
		R1928 = mro_crism_lookupwv(1928,wvt)
		R1934 = mro_crism_lookupwv(1934,wvt)
		R1941 = mro_crism_lookupwv(1941,wvt)
		R1862 = mro_crism_lookupwv(1862,wvt)
		R1869 = mro_crism_lookupwv(1869,wvt)
		R1875 = mro_crism_lookupwv(1875,wvt)
		R2112 = mro_crism_lookupwv(2112,wvt)
		R2120 = mro_crism_lookupwv(2120,wvt)
		R2126 = mro_crism_lookupwv(2126,wvt)
		pcube(*,*,paramidx)= 1.0-(((cube(*,*,R1908)+cube(*,*,R1914)+cube(*,*,R1921)+cube(*,*,R1928)+cube(*,*,R1934)+cube(*,*,R1941))) / $
						(cube(*,*,R1862) + cube(*,*,R1869)+cube(*,*,R1875)+cube(*,*,R2112)+cube(*,*,R2120)+cube(*,*,R2126) ) )
	endif


	;------------------------------------------------------------------------
	; Find the 1.98 micron sulfite band depth (BD1980): (wvs 1980-1987 (in), 1921, 2040)
	;------------------------------------------------------------------------
	paramidx++
	if (index_vec[paramidx] EQ 1) then begin
		R1921 = mro_crism_lookupwv(1921,wvt)
		R1980 = mro_crism_lookupwv(1980,wvt)
		R1987 = mro_crism_lookupwv(1987,wvt)
		R2040 = mro_crism_lookupwv(2040,wvt)
		pcube[*,*,paramidx] = $
			1.0-(((cube[*,*,R1980]+cube[*,*,R1987])) / (cube[*,*,R1921] + cube[*,*,R2040] ) )
	endif


	;------------------------------------------------------------------------
	;Find the 2.20 micron AL-OH band depth (BD2200): (wvs 2132, 2146, 2199, 2205, 2252, 2258)
	;------------------------------------------------------------------------
	paramidx++
	if (index_vec[paramidx] EQ 1) then begin
		R2132 = mro_crism_lookupwv(2132,wvt)
		R2146 = mro_crism_lookupwv(2146,wvt)
		R2199 = mro_crism_lookupwv(2199,wvt)
		R2205 = mro_crism_lookupwv(2205,wvt)
		R2252 = mro_crism_lookupwv(2252,wvt)
		R2258 = mro_crism_lookupwv(2258,wvt)
		pcube[*,*,paramidx] = 1.0-(((cube(*,*,R2199)+cube(*,*,R2205))*2.0) / $
                            (cube(*,*,R2132) + cube(*,*,R2146)+cube(*,*,R2252)+cube(*,*,R2258) ) )
	endif

	;------------------------------------------------------------------------
	; Find the 2.22 and 2.28 micron doublet band depth
	; (Doub2200): (wvs 2172, 2205, 2258, 2311)
	;------------------------------------------------------------------------
	paramidx++
	if (index_vec[paramidx] EQ 1) then begin
		R2172 = mro_crism_lookupwv(2172,wvt)
		R2205 = mro_crism_lookupwv(2205,wvt)
		R2258 = mro_crism_lookupwv(2258,wvt)
		R2311 = mro_crism_lookupwv(2311,wvt)
		pcube[*,*,paramidx] = 1.0-((cube(*,*,R2205)+cube(*,*,R2258)) / (cube(*,*,R2172)+cube(*,*,R2311)) )
	endif

	;------------------------------------------------------------------------
	; Find the 2.23 micron band depth  (K. Lichtenberg)
	;------------------------------------------------------------------------
	paramidx++
	if (index_vec[paramidx] EQ 1) then begin
		R2231 = mro_crism_lookupwv(2231,wvt)                            ; 
		R2258 = mro_crism_lookupwv(2258,wvt)                            ; 
		R2251 = mro_crism_lookupwv(2251,wvt) 
		R2212 = mro_crism_lookupwv(2212,wvt)
		R2198 = mro_crism_lookupwv(2198,wvt)
 
		WL1 = mro_crism_lookupwv(2212,wvt,/w)
		WL2 = mro_crism_lookupwv(2198,wvt,/w) 
		WC = mro_crism_lookupwv(2231,wvt,/w)
		WH1 = mro_crism_lookupwv(2258,wvt,/w)
		WH2 = mro_crism_lookupwv(2251,wvt,/w)

		b = (WC(0)-(WL1(0)+WL2(0))/2)/((WH1(0)+WH2(0))/2-(WL1(0)+WL2(0))/2)   ; 
		a = 1-b  

		pcube(*,*,paramidx) = 1-( cube(*,*,R2231) / ( a*(cube(*,*,R2212)+ cube(*,*,R2198))/2 + b *(cube(*,*,R2258)+ cube(*,*,R2251))/2 ) )
	endif
 
	;------------------------------------------------------------------------
	;Find the 2.50 micron band depth (BD2500): (wvs 2380, 2500, 2510, 2540)
	;------------------------------------------------------------------------
	paramidx++
	if (index_vec[paramidx] EQ 1) then begin
		R2380 = mro_crism_lookupwv(2380,wvt)
		R2500 = mro_crism_lookupwv(2500,wvt)
		R2510 = mro_crism_lookupwv(2510,wvt)
		R2540 = mro_crism_lookupwv(2540,wvt)
		pcube(*,*,paramidx) = 1.0-( (cube(*,*,R2500) + cube(*,*,R2510)) / ( cube(*,*,R2540) + cube(*,*,R2380) ) )
	endif


   if (det eq 2) then pcubej(*,*,NUM_VNIR_PARAM:*)=pcube

endif

if (det eq 2) then pcube=pcubej

return, pcube[*,*,where(index_vec_in EQ 1)]

end
