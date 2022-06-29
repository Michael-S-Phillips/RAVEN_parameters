;ABSTRACT
;Function to calculate summary products for CRISM. 
;
;INPUTS
;   cube = input cube -> dark subtracted, I/F cube or proxy.
;           X=along track (640),y = number of frames, z= #wavelengths
;           so, a given wvln slice of the cube is (*,*,wvln)
;   det  = detector flag. default = VNIR
;   wvlut = filename of wavelength lookup table. Default is
;   wvt = wavelength array passed from ENVI header

;HISTORY
; V1.0 Created 2005/05/10 by Noam Izenberg (NRI), JHU/APL
; V2.0 Began incorporating OMEGA solutions/subroutines from Shannon Pelkey, 2005/09/01 NRI
;		- csumprod_rpeak1 function
; V3.0 Replace product with pcube do all VNIR or all IR slices at once, and return a product cube.
; v3.1  2006/04/07 Updating to Zeta2 wavelength table and Shannon's version 8 getparams
;	Replaced static wavelength with lookupwv.pro
;	Finishing the remaining curve fit solutions
;Notes
; lookupwv inputs are: target wvln in nm, detector (0=VNir, 1=IR, and lut variable
;** Wavelength Table Zeta2 summary (row number to wavelength for both detectors)
;	 is included in Appendix after program end.
;** Channel math spreadsheet is used to determine a,b,c,d,e weighting parameters
;** In order to Match compatability with OMEGA parameters, csumprod uses wavelengths as
; close to those in the CRISM SIS and cal report (10) as possible where They do not
; conflict with "bad" CRSIM rows (i.e. covered rows or rows with
; filter seams).
;
; Jun 27, 2007: S.Pelkey
; --modified formulations to the following parameters:
; VNIR:  R410 became R440, BD530, SH600, BD640, BD860, BD920
; IR: VAR, BD1500, BD1750, BD1900, BDI2000,OLINDEX,LCPINDEX,HCPINDEX,IRR3 
; --modified to process joined detectors (like TILES)
;-----------------------------------------------------------------
;
;10/05/2007 (fps)
;   Forked development of summary parameter function
;   Renamed to mro_crism_summary_params_dev to avoid conflicts with released version
;   Added igore_val keyword, ignore_val masking of return array
;10/08/2007 (fps)
;   Added option to return summary parameter band names instead of summary parameter cube (ad-hoc keyword override)
;03/04/2008 (fps)
;   Total overhaul of RPEAK1, BDI1000VIS
;   Modified SH600
;   Added informative console output to the modified calculations
;03/07/2008 (fps)
;   Fixed ISLOPE1 param name in the band_names vectors
;03/19/2008 (fps)
;   Fixed coding error in BDI1000VIS to integrate (1.0 - normalized values) rather than (1.0 - values)
;   Total overhaul of BDI1000IR
;03/23/2008 (fps)
;   Total overhaul of BDI2000IR
;03/31/2008  (fps)
;  Significant RPEAK1 speed enhancement - use IDL poly function instead of custom polynomial calculation loop - remember the temporary directive!
; 05/06/2008 (mfm)
;   Modified D2300 and BDCARB to use 2120 instead of 2140. This is per recomendation of F. Seelos 
;   and B. Ehlmann.
; 03/10/2009 (mfm)
;	Use cat_int_tabulated() in place of int_tabulated() to avoid crash on repeated x values
;	in some cases when bad bands specified.
; 10/25/2010 (mfm)
;	Fix error in BD2290 - continuum interpolation weights were hardcoded and reversed
; 10/26/2010 (mfm)
;	Fix 3 errors in BD2000co2:
;	    1) Wavelength for WH is LAMB=2170, not 3390
;	    2) continuum interpolation weights were hardcoded and reversed
;	    3) a weights cube(WH=2170), b weights cube(WL=1815); not the other way around
; 11/04/2010 (fps)
;   Minor udpate to CRISM NAN handling to allow for input cubes that do not have any CRISM NANs - e.g. Dick and Sandra's joined DISORT ratio cubes
; 12/20/2010 (mfm)
;   Resolved svn conflict with last updates (mfm vs fps) - no real problem

function mro_crism_summary_params, cube, det, wvt, index_vec=index_vec, ignore_val = ignore_val, band_names = band_names

if (not(keyword_set(ignore_val))) then begin
    ignore_val = 65535.0
endif

if (not(keyword_set(index_vec))) then begin
    if (det EQ 0) then index_vec = replicate(1, 11) ;S
    if (det EQ 1) then index_vec = replicate(1, 34) ;L
    if (det EQ 2) then index_vec = replicate(1, 45) ;J
endif

if (keyword_set(band_names)) then begin

    case(det) of
    	0:begin
	    band_names = ['R770', 'RBR', 'BD530', 'SH600', 'BD640', 'BD860', 'BD920', 'RPEAK1', 'BDI1000VIS', $
	    	    	  'R440', 'IRR1']
	end
	
	1:begin
	    band_names = ['BD11000IR', 'IRA', 'OLINDEX', 'LCPINDEX', 'HCPINDEX', 'VAR', 'ISLOPE1', 'BD1435', 'BD1500', 'ICER1', $
	    	    	  'BD1750', 'BD1900', 'BDI2000', 'BD2100', 'BD2210', 'BD2290', 'D2300', 'SINDEX', 'ICER2', 'BDCARB', $
			  'BD3000', 'BD3100', 'BD3200', 'BD3400', 'CINDEX', $
			  'BD1270O2', 'BD1400H2O', 'BD2000CO2', 'BD2350', 'BD2600', 'IRR2', 'R2700', 'BD2700', 'IRR3']
	end
	
	2:begin
	    band_names = ['R770', 'RBR', 'BD530', 'SH600', 'BD640', 'BD860', 'BD920', 'RPEAK1', 'BDI1000VIS', $
	    	    	  'R440', 'IRR1', $
			  'BD11000IR', 'IRA', 'OLINDEX', 'LCPINDEX', 'HCPINDEX', 'VAR', 'ISLOPE1', 'BD1435', 'BD1500', 'ICER1', $
	    	    	  'BD1750', 'BD1900', 'BDI2000', 'BD2100', 'BD2210', 'BD2290', 'D2300', 'SINDEX', 'ICER2', 'BDCARB', $
			  'BD3000', 'BD3100', 'BD3200', 'BD3400', 'CINDEX', $
			  'BD1270O2', 'BD1400H2O', 'BD2000CO2', 'BD2350', 'BD2600', 'IRR2', 'R2700', 'BD2700', 'IRR3']
	end
    
    	else:begin
    	   band_names =  -1
	endelse
    endcase

    
    return_band_names = band_names[where(index_vec EQ 1)]

    return, return_band_names

endif


sz = size(cube)                  ; x = spatial (640 detector columns, unless binned)
		 		 ; y = spatial (# of frames along track)
		 		 ; z = spectral # wavelengths 

NB = sz[3]       ; spectral bands
NX = sz[1]	 ; spatial detector
NY = sz[2]	 ; spatial along track


;Calculate ignore value mask band
cube_median = median(cube, dimension = 3)

cube_median_hist = histogram(cube_median, min = ignore_val, max = ignore_val, reverse_indices = ri_cube_median_hist)
cube_mask_band = make_array(nx, ny, value = 0.0)
if (cube_median_hist[0] GT 0) then begin
    cube_mask_band[ri_cube_median_hist[ri_cube_median_hist[0]:ri_cube_median_hist[1]-1]] = ignore_val
endif


;dname=['VNIR','IR','Joined']
;Detector = dname[det]


; For tiles:
if (det eq 2) then pcubej = make_array(nx, ny, 45, value = 65535.)

; Make a copy of the input index vector:
index_vec_in=index_vec

;stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PARAMETERS USING VNIR WAVELENGTHS:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (det eq 0) or (det eq 2) then begin  
;print, 'VNIR Parameter Calculation...'
;print, systime(0)
vnir_start_time = systime(1)

        ; Set index_vec just to VNIR params:
        if (det eq 2) then index_vec=index_vec_in[0:10]

        pcube = make_array(nx, ny, 11, value = 65535.)

	;------------------------------------------------------------------------
	; Find the 0.77 micron reflectance (R770):
	;------------------------------------------------------------------------
	if (index_vec[0] EQ 1) then begin
		R770 = mro_crism_lookupwv(770,wvt)
;                w=where(cube(*,*,R770) ne 65535.)
;	 	if (w[0] ne -1) then pcube(w,*,0) = (cube(*,*,R770))[w]
		pcube(*,*,0) = cube(*,*,R770)
	endif

	;------------------------------------------------------------------------
	; Find the red/blue ratio (RBR):
	;------------------------------------------------------------------------
	if (index_vec[1] EQ 1) then begin
		if (n_elements(R770) EQ 0) then begin
			R770 = mro_crism_lookupwv(770,wvt)
		endif

		R440 = mro_crism_lookupwv(440,wvt)
		pcube(*,*,1) = cube(*,*,R770) / cube(*,*,R440)
	endif

	;------------------------------------------------------------------------
	; Find the 0.53 micron band depth (BD530):
	; SLM modified 26 Jun 2007 to move long wavelength shoulder out 
	; of VNIR filter zone boundary, from 648 to 709 nm
	;------------------------------------------------------------------------
	if (index_vec[2] EQ 1) then begin
		if (n_elements(R440) EQ 0) then begin
			R440 = mro_crism_lookupwv(440,wvt)
		endif
		R530 = mro_crism_lookupwv(530,wvt)
		R709 = mro_crism_lookupwv(709,wvt)
		WL = mro_crism_lookupwv(440,wvt,/w)
		WC = mro_crism_lookupwv(530,wvt,/w)
		WH = mro_crism_lookupwv(709,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))		; a gets multipled by the longer band
		b = 1-a		      			; b gets multiplied by the shorter band
		pcube(*,*,2) = 1 - ( cube(*,*,R530) / (a * cube(*,*,R709) + b * cube(*,*,R440) ) )
            endif

	;------------------------------------------------------------------------
	; Find the 0.60 micron shoulder height (SH600):
	; SLM modified 26 Jun 2007 to move long wavelength shoulder out 
	; of VNIR filter zone boundary, from 680 to 709 nm
	; 03/04/2008	(fps)
	;   Changed shoulder edge wavelengths reference to [533,710] from [530,709] to resolve the same detector rows in multi- and hyper-spectral data
	;   Reformulated as an 'inverted' band depth 
	;------------------------------------------------------------------------
	if (index_vec[3] EQ 1) then begin
        if (n_elements(R533) EQ 0) then begin
            R533 = mro_crism_lookupwv(533,wvt)
        endif
	    R600 = mro_crism_lookupwv(600,wvt)
	    R710 = mro_crism_lookupwv(710,wvt)
	    WL = mro_crism_lookupwv(533,wvt,/w)
	    WC = mro_crism_lookupwv(600,wvt,/w)
	    WH = mro_crism_lookupwv(710,wvt,/w)
	    a = (WC(0)-WL(0))/(WH(0)-WL(0)) 	; a gets multipled by the longer band
	    b = 1.0-a		      	    	; b gets multiplied by the shorter band
	    ;pcube(*,*,3) = cube(*,*,R600) / (b * cube(*,*,R533) + a * cube(*,*,R710))
	    pcube[*,*,3] = 1.0 - ((b * cube[*,*,R533] + a * cube[*,*,R710]) / cube[*,*,R600])
    	    
	    ;print, systime(0)
	endif

	;------------------------------------------------------------------------
	; Find the 0.64 micron band depth (BD640): (really calculated at 0.648)
	; SLM modified 26 Jun 2007 to move long wavelength shoulder out 
	; of VNIR filter zone boundary, from 680 to 709 nm
	; NOTE: THIS PARAMETER WILL BE DEGRADED BECAUSE THE KEY WAVELENGTH NEAR
	; 648 NM IS LOCATION WITHIN THE VNIR FILTER ZONE BOUNDARY
	;------------------------------------------------------------------------
	if (index_vec[4] EQ 1) then begin
		if (n_elements(R600) EQ 0) then begin
			R600 = mro_crism_lookupwv(600,wvt)
		endif
		if (n_elements(R709) EQ 0) then begin
			R709 = mro_crism_lookupwv(709,wvt)
		endif
		if (n_elements(R648) EQ 0) then begin
			R648 = mro_crism_lookupwv(648,wvt)
		endif
		WL = mro_crism_lookupwv(600,wvt,/w)
		WC = mro_crism_lookupwv(648,wvt,/w)
		WH = mro_crism_lookupwv(709,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer band
		b = 1-a		      		  ; b gets multiplied by the shorter band
		pcube(*,*,4) = 1 - ( cube(*,*,R648) /  (b * cube(*,*,R600) + a * cube(*,*,R709) ) )
        endif

	;------------------------------------------------------------------------
	; Find the 0.86 micron band depth (BD860): ('hematite band')
	; SLM modified 26 Jun 2007 to move long wavelength shoulder 
	; from 920 to 984 nm to sample more spectral curvature
	;------------------------------------------------------------------------
	if (index_vec[5] EQ 1) then begin
		R800 = mro_crism_lookupwv(800,wvt)
		R860 = mro_crism_lookupwv(860,wvt)
		R984 = mro_crism_lookupwv(984,wvt)
		WL = mro_crism_lookupwv(800,wvt,/w)
		WC = mro_crism_lookupwv(860,wvt,/w)
		WH = mro_crism_lookupwv(984,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer band
		b = 1-a		      				  ; b gets multiplied by the shorter band
		pcube(*,*,5) = 1 - ( cube(*,*,R860) /  (b * cube(*,*,R800) + a * cube(*,*,R984) ) )
        endif

	;------------------------------------------------------------------------
	; Find the 0.92 micron band depth (BD920): ('Pseudo BDI1000 VIS')
	; SLM modified 26 Jun 2007 to move long wavelength shoulder 
	; from 1020 to 984 nm to avoid calibration inaccuracy at VNIR wavelengths
	; greater than 1010 nm
	;------------------------------------------------------------------------
	if (index_vec[6] EQ 1) then begin
		R800 = mro_crism_lookupwv(800,wvt)
		R920 = mro_crism_lookupwv(920,wvt)
		R984 = mro_crism_lookupwv(984,wvt)
		WL = mro_crism_lookupwv(800,wvt,/w)
		WC = mro_crism_lookupwv(920,wvt,/w)
		WH = mro_crism_lookupwv(984,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))	; a gets multipled by the longer band
		b = 1-a		      				; b gets multiplied by the shorter band
		pcube(*,*,6) = 1 - ( cube(*,*,R920) /  (b * cube(*,*,R800) + a * cube(*,*,R984) ) )
	endif

	;------------------------------------------------------------------------
	; Find reflectance peak 1 (RPEAK1):
    	; 03/04/2008	(fps) 
	;   Total rewrite of the RPEAK1 calculation
	;   Simple poyynomial fit of reflectances over specified wavelength range
	;   Wavelength vector ~ [442,533,600,710,740,775,800,833,860,892] nm - Selected to resolve to same detector rows in both multi- and hyper-spectral cubes
	;   Bypassing ~ 660nm filter boundary as much as possible
	;   Polynomial order reduced from 5 to 4 - compromise between bining of results (e.g. spectral smile effect) and representation of spectral shape
	;   Results are much improved, but RPEAK1 run time is ~4 minutes for an FRT on crism-lnx## machine
	; 03/31/2008  (fps)
	;   Dramatic speed improvement using IDL poly function (improved loop data handling)
	;------------------------------------------------------------------------
	if ((index_vec[7] EQ 1) OR (index_vec[8] EQ 1)) then begin
    	    ;print, 'RPEAK1'
	    ;print, systime(0)
	    
    	    show_plot = 0

    	    if (n_elements(R442) EQ 0) then begin
		R442 = mro_crism_lookupwv(442,wvt)
	    endif
    	    if (n_elements(R533) EQ 0) then begin
		R533 = mro_crism_lookupwv(533,wvt)
	    endif
    	    if (n_elements(R600) EQ 0) then begin
		R600 = mro_crism_lookupwv(600,wvt)
	    endif
	    if (n_elements(R710) EQ 0) then begin
		R710 = mro_crism_lookupwv(710,wvt)
	    endif
	    if (n_elements(R740) EQ 0) then begin
		R740 = mro_crism_lookupwv(740,wvt)
	    endif
	    if (n_elements(R775) EQ 0) then begin
		R775 = mro_crism_lookupwv(775,wvt)
	    endif
	    if (n_elements(R800) EQ 0) then begin
		R800 = mro_crism_lookupwv(800,wvt)
	    endif
	    if (n_elements(R833) EQ 0) then begin
		R833 = mro_crism_lookupwv(833,wvt)
	    endif
	    if (n_elements(R860) EQ 0) then begin
		R860 = mro_crism_lookupwv(860,wvt)
	    endif
	    if (n_elements(R892) EQ 0) then begin
		R892 = mro_crism_lookupwv(892,wvt)
	    endif
	    if (n_elements(R925) EQ 0) then begin
		R925 = mro_crism_lookupwv(925,wvt)
	    endif

    	    wavelength_indx = [R442,R533,R600,R710,R740,R775,R800,R833,R860,R892,R925]
    	    wavelength_vec = wvt[wavelength_indx]

    	    num_model_points = 1001
	    poly_degree = 4.0
    	    model_wavelength_vec = findgen(num_model_points) / (num_model_points - 1.0) * (max(wavelength_vec) - min(wavelength_vec)) + min(wavelength_vec)

    	    ;print, '  Wavelength vector:'
    	    ;print, wavelength_vec
	    ;print, '  Index vector:'
    	    ;print, wavelength_indx

    	    ;Define results band
	    rpeak_wavelength = make_array(nx,ny, value = ignore_val, /float)
	    rpeak_value = make_array(nx,ny, value = ignore_val, /float)

    	    ;Isolate appropriate bands from input_cube

    	    rpeak_cube = cube[*,*,wavelength_indx]
	    
	    if (show_plot EQ 1) then begin
            	window, 17, xsize = 600, ysize = 600, retain = 2
	    endif
	    
    	    for k = 0, nx -1 do begin
	    	;print, k, systime(0)
	    	for l = 0, ny -1 do begin
		    spec_vec = rpeak_cube[k,l,*]
		    
		    if (total(spec_vec EQ ignore_val) EQ 0) then begin
    	    	    	
			if (show_plot EQ 1) then begin
		    	    plot, wavelength_vec, spec_vec, psym = -1, $
			      	yrange = [0.0, max(spec_vec) + 0.10 * max(spec_vec)], ysty=1, $
			      	xtitle = 'Wavelength (nm)', ytitle = 'Spectral value', $
			      	title = 'Pixel x=' + strtrim(string(k),2) + ' y=' + strtrim(string(l),2)	
			endif
			
			rpeak_params = poly_fit(wavelength_vec, spec_vec, poly_degree)
			
			;model_spec = make_array(num_model_points, value = 0.0)
			;for m = 0, poly_degree do begin
			;    model_spec = model_spec + rpeak_params[m] * model_wavelength_vec^(float(m))
			;endfor
			model_spec = poly(model_wavelength_vec, rpeak_params)
			
			if (show_plot EQ 1) then begin
			    oplot, model_wavelength_vec, model_spec
       	    	    	    ;wait, 0.25
			endif

    	    	    	model_spec_max = max(model_spec, model_spec_max_indx)
    	    	    	model_spec_max_wavelength = model_wavelength_vec[model_spec_max_indx]

			rpeak_wavelength[k,l] = model_spec_max_wavelength
			rpeak_value[k,l] = model_spec_max
			
		    endif
		endfor
	    endfor

	pcube[*,*,7] = rpeak_wavelength / 1000.0
	
	;print, systime(0)
	
	endif


	;------------------------------------------------------------------------
	; Find the 1 micron band depth (BDI1000VIS):
	; (wait on 1 micron band center) VIS and IR
	; 03/04/2008	(fps)
	;   Total rewrite of BDI1000VIS calculation to take advantage of revised RPEAK1 procedure
	;   Wavelength vector ~ [833,860,892,925,951,984,1023] nm - Selected to resolve to same detector rows in both multi- and hyper-spectral cubes
	;   Integration is performed in microns to keep results in range of simple BD parameters [0,1]
	;------------------------------------------------------------------------
	if (index_vec[8] EQ 1) then begin
    	    ;print, 'BDI1000VIS'
	    ;print, systime(0)
	    
    	    if (n_elements(R833) EQ 0) then begin
    	    	R833 = mro_crism_lookupwv(833,wvt)
    	    endif
    	    if (n_elements(R860) EQ 0) then begin
    	    	R860 = mro_crism_lookupwv(860,wvt)
    	    endif
    	    if (n_elements(R892) EQ 0) then begin
	    	R892 = mro_crism_lookupwv(892,wvt)
	    endif
	    if (n_elements(R925) EQ 0) then begin
    	    	R925 = mro_crism_lookupwv(925,wvt)
	    endif
	    if (n_elements(R951) EQ 0) then begin
    	    	R951 = mro_crism_lookupwv(951,wvt)
    	    endif
    	    if (n_elements(R984) EQ 0) then begin
	    	R984 = mro_crism_lookupwv(984,wvt)
	    endif
	    if (n_elements(R1023) EQ 0) then begin
    	    	R1023 = mro_crism_lookupwv(1023,wvt)
	    endif
    	    
    	    wavelength_indx = [R833,R860,R892,R925,R951,R984,R1023]	
      	    wavelength_vec = wvt[wavelength_indx]
	    
	    wavelength_vec_um = wavelength_vec / 1000.0
	    
	    ;print, '  Wavelength vector:'
	    ;print, wavelength_vec
	    ;print, '  Index vector:'
    	    ;print, wavelength_indx		

    	    bdi1000vis_value = make_array(nx, ny, value = ignore_val)
	    
    	    bdi1000_cube = cube[*,*,wavelength_indx]
    	    rpeak_value_cube = make_array(size = size(bdi1000_cube), value = ignore_val)
    	    for j = 0, n_elements(wavelength_indx) -1 do begin
	    	rpeak_value_cube[*,*,j] = rpeak_value
	    endfor

    	    bdi1000_normalized_cube = bdi1000_cube / rpeak_value_cube


    	    for k = 0, nx -1 do begin
	    	for l = 0, ny -1 do begin

    	    	    spec_vec = bdi1000_normalized_cube[k,l,*]
    	    	    check_vec = bdi1000_cube[k,l,*]
		    
		    if (total(check_vec EQ ignore_val) EQ 0) then begin
		    
		    	bdi1000vis_value[k,l] = cat_int_tabulated(wavelength_vec_um, 1.0 - spec_vec)
		
		    endif
		    
		endfor
	    endfor

    	    pcube[*,*,8] = bdi1000vis_value

	    ;print, systime(0)
	
	endif

	;************************************************************************
	;PARAMETERS FOR VNIR NON-ATMOSPHERICALLY-CORRECTED CALIBRATED DATA
	;************************************************************************

	;------------------------------------------------------------------------
	;Find the 0.44 micron reflectance R440 (changed from R410):
	; SLM modified 26 Jun 2007 to move from 410 to 440 nm 
	; to avoid calibration errors from VNIR scattered light removal -
	; the lowest 'valid' wavelength varies from 410 nm in typical terrain
	; are about 440 nm in polar ice regions
	;------------------------------------------------------------------------
	if (index_vec[9] EQ 1) then begin
		R440 = mro_crism_lookupwv(440,wvt)
		pcube(*,*,9) = cube(*,*,R440)
	endif

	;------------------------------------------------------------------------
	;Find the IR ratio 1 (IRR1):
	;------------------------------------------------------------------------
	if (index_vec[10] EQ 1) then begin
		if (n_elements(R800) EQ 0) then begin
			R800 = mro_crism_lookupwv(800,wvt)
		endif

		R1020 = mro_crism_lookupwv(1020,wvt)
		pcube(*,*,10) = cube(*,*,R800) / cube(*,*,R1020)
        endif


    if (det eq 2) then pcubej(*,*,0:10)=pcube
    
endif                ; VNIR or J case

        

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PARAMETERS USING IR WAVELENGTHS:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (det eq 1) or (det eq 2) then begin

    ;print, 'IR Parameter Calculation...'
    ;print, systime(0)
    ir_start_time = systime(1)
    
        ; Reset index_vec to IR params:
        if (det eq 2) then index_vec=index_vec_in[11:*]

	pcube = make_array(nx, ny, 34, value = 65535.)



	;------------------------------------------------------------------------
	; Find the 1 micron band depth (BDI1000IR):
	; 03/19/2008	(fps) 
	;   Total rewrite of the BDI1000IR calculation
	;   Revised continuum calculation - used for both BDI1000IR and BDI2000
	;   	Find median value at multispectral wavelengths over the interval ~[1330,1810]
	;   	Find median value at multispectral wavelengths over the interval ~[2430,2600]
	;   	Linear fit between the median (wavelength, value) points is the continuum
    	;   Revised ~1000 nm integrated band depth calculation
	;   	Integration of normalized reflectance values over the inteveral ~[1045,1255]
	; 04/02/2008	(fps)
	;   Refined continuum calculation used in BDI1000IR and BDI2000 calculations
	;   	Find upper quartile (75th percentile) value at multispectral wavelengths over the interval ~[1330,1810]
	;   	Linear fit between the ~2500 nm median (wavelength, value) and ~1500 nm upper quartile (wavelength, value) points is the continuum
	;   	  
	;------------------------------------------------------------------------
	if ((index_vec[0] EQ 1) OR (index_vec[12] EQ 1)) then begin
    	
    	;print, 'BDI1000IR'
    	;print, systime(0)

    	show_map = 0
    	show_plot = 0
	
    	;Contunuum calculation used in BDI1000IR and BDI2000
	    
    	    if (n_elements(R1275) EQ 0) then begin
    	    	R1275 = mro_crism_lookupwv(1275,wvt)
    	    endif
    	    if (n_elements(R1330) EQ 0) then begin
    	    	R1330 = mro_crism_lookupwv(1330,wvt)
    	    endif
    	    if (n_elements(R1370) EQ 0) then begin
	    	R1370 = mro_crism_lookupwv(1370,wvt)
	    endif
	    if (n_elements(R1395) EQ 0) then begin
    	    	R1395 = mro_crism_lookupwv(1395,wvt)
	    endif
	    if (n_elements(R1425) EQ 0) then begin
    	    	R1425 = mro_crism_lookupwv(1425,wvt)
    	    endif
    	    if (n_elements(R1465) EQ 0) then begin
	    	R1465 = mro_crism_lookupwv(1465,wvt)
	    endif
	    if (n_elements(R1500) EQ 0) then begin
    	    	R1500 = mro_crism_lookupwv(1500,wvt)
	    endif
	    if (n_elements(R1505) EQ 0) then begin
    	    	R1505 = mro_crism_lookupwv(1505,wvt)
	    endif
	    if (n_elements(R1560) EQ 0) then begin
    	    	R1560 = mro_crism_lookupwv(1560,wvt)
	    endif
	    if (n_elements(R1625) EQ 0) then begin
    	    	R1625 = mro_crism_lookupwv(1625,wvt)
	    endif
	    if (n_elements(R1660) EQ 0) then begin
    	    	R1660 = mro_crism_lookupwv(1660,wvt)
	    endif
	    if (n_elements(R1690) EQ 0) then begin
    	    	R1690 = mro_crism_lookupwv(1690,wvt)
	    endif
	    if (n_elements(R1750) EQ 0) then begin
    	    	R1750 = mro_crism_lookupwv(1750,wvt)
	    endif
	    if (n_elements(R1810) EQ 0) then begin
    	    	R1810 = mro_crism_lookupwv(1810,wvt)
	    endif
	    if (n_elements(R1875) EQ 0) then begin
    	    	R1875 = mro_crism_lookupwv(1875,wvt)
	    endif

;    	    wavelength_indx = [R1275,R1330,R1370,R1395,R1425,R1465,R1500,R1505,R1560,R1625,R1660,R1690,R1750,R1810,R1875]
;	    wavelength_indx = [R1330,R1370,R1395,R1425,R1465,R1500,R1505,R1560,R1625,R1660,R1690,R1750,R1810,R1875]
	    wavelength_indx = [R1330,R1370,R1395,R1425,R1465,R1500,R1505,R1560,R1625,R1660,R1690,R1750,R1810]
    	    wavelength_vec = wvt[wavelength_indx]
	    
    	    ;print, '  Wavelength vector:'
    	    ;print, wavelength_vec
	    ;print, '  Index vector:'
    	;    print, wavelength_indx
	    
	    cont1_cube = cube[*,*,wavelength_indx]
	    ;help, cont1_cube

    	    percentile1_frame = make_array(nx, ny, value = ignore_val)
	    percentile1_wavelength = make_array(nx,ny, value = ignore_val)
    	    ;help, percentile1_frame    	    
	    ;help, percentile1_wavelength

    	    for i = 0, nx -1 do begin
	    	for j = 0, ny -1 do begin
		
		    spec_vec = cont1_cube[i,j,*]
		    
		    if (total(spec_vec EQ ignore_val) EQ 0) then begin
		    
		    	sort_indx = sort(spec_vec)
		    	select_indx = sort_indx[n_elements(wavelength_indx) * 0.75] ;~upper quartile indx
		    
		    	percentile1_frame[i,j] = spec_vec[select_indx]
		    	percentile1_wavelength[i,j] = wavelength_vec[select_indx]
			
		    endif
		    		
		endfor
	    endfor
	    	    
;	    median1_frame = median(cont1_cube, dimension = 3)
;	    help, median1_frame
;	    median1_wavelength = median(wavelength_vec, /even)
;	    help, median1_wavelength
	    	    
	    if (n_elements(R2430) EQ 0) then begin
    	    	R2430 = mro_crism_lookupwv(2430,wvt)
    	    endif
    	    if (n_elements(R2460) EQ 0) then begin
    	    	R2460 = mro_crism_lookupwv(2460,wvt)
    	    endif
    	    if (n_elements(R2530) EQ 0) then begin
	    	R2530 = mro_crism_lookupwv(2530,wvt)
	    endif
    	    if (n_elements(R2600) EQ 0) then begin
	    	R2600 = mro_crism_lookupwv(2600,wvt)
	    endif
	    
	    wavelength_indx = [R2430,R2460,R2530,R2600]
	    wavelength_vec = wvt[wavelength_indx]
	    	    	    
    	    ;print, '  Wavelength vector:'
    	    ;print, wavelength_vec
	    ;print, '  Index vector:'
    	 ;   print, wavelength_indx
	    
    	    cont2_cube = cube[*,*,wavelength_indx]
	    ;help, cont2_cube
	    
	    median2_frame = median(cont2_cube, dimension = 3)
	    ;help, median2_frame
	    
	    median2_wavelength = median(wavelength_vec, /even)
	    ;help, median2_wavelength
	    
    	    cont_slope_frame = (median2_frame - percentile1_frame) / (median2_wavelength - percentile1_wavelength)

    	    if (keyword_set(show_map)) then begin	    
		window, 10, xsize = 700, ysize = 500, retain = 2, title = 'Percentile1 Frame'
	    	tv, bytscl(percentile1_frame, min = 0.05, max = 0.35)
;        	window, 10, xsize = 700, ysize = 500, retain = 2, title = 'Median1 Frame'
;	    	tv, bytscl(median1_frame, min = 0.05, max = 0.35)

    	    	window, 11, xsize = 700, ysize = 500, retain = 2, title = 'Percentile1 Wavelength'
	    	tv, bytscl(percentile1_wavelength, min = 1300, max = 1800)
		
    	    	window, 12, xsize = 700, ysize = 500, retain = 2, title = 'Median2 Frame'
	    	tv, bytscl(median2_frame, min = 0.05, max = 0.35)

    	    	window, 13, xsize = 700, ysize = 500, retain = 2, title = 'Continuum Slope'
	    	tv, bytscl(cont_slope_frame, min = min(cont_slope_frame), max = max(cont_slope_frame))
    	    endif

    	    ;Normalize and integrate ~1 um band

    	    if (n_elements(R1020) EQ 0) then begin
    	    	R1020 = mro_crism_lookupwv(1020,wvt)
    	    endif
    	    if (n_elements(R1045) EQ 0) then begin
    	    	R1045 = mro_crism_lookupwv(1045,wvt)
    	    endif
    	    if (n_elements(R1080) EQ 0) then begin
	    	R1080 = mro_crism_lookupwv(1080,wvt)
	    endif
	    if (n_elements(R1150) EQ 0) then begin
    	    	R1150 = mro_crism_lookupwv(1150,wvt)
	    endif
	    if (n_elements(R1210) EQ 0) then begin
    	    	R1210 = mro_crism_lookupwv(1210,wvt)
	    endif
	    if (n_elements(R1250) EQ 0) then begin
    	    	R1250 = mro_crism_lookupwv(1250,wvt)
	    endif
	    if (n_elements(R1255) EQ 0) then begin
    	    	R1255 = mro_crism_lookupwv(1255,wvt)
	    endif
;    	    wavelength_indx = [R1020,R1045,R1150,R1210]
	    wavelength_indx = [R1045,R1150,R1210,R1250,R1255]
	    wavelength_vec = wvt[wavelength_indx]

    	    wavelength_vec_um = wavelength_vec / 1000.0

	    ;print, '  Wavelength vector:'
    	 ;   print, wavelength_vec
	   ; print, '  Index vector:'
    	;    print, wavelength_indx

    	    bdi1000ir_cube = cube[*,*,wavelength_indx]
    	    
	    cont_cube = make_array(nx, ny, n_elements(wavelength_vec), value = ignore_val)

	    for k = 0, n_elements(wavelength_vec) -1 do begin
    	    	cont_cube[*,*,k] = cont_slope_frame * (wavelength_vec[k] - median2_wavelength) + median2_frame
	    endfor

    	    bdi1000ir_normalized_cube = bdi1000ir_cube / cont_cube

    	    bdi1000ir_value = make_array(nx, ny, value = ignore_val)

    	    if (keyword_set(show_plot)) then begin
	        window, 0, xsize = 1000, ysize = 500, retain=2
	    endif

    	    for k = 0, nx -1 do begin
	    	for l = 0, ny -1 do begin

    	    	    spec_vec = bdi1000ir_normalized_cube[k,l,*]
    	    	    check_vec = bdi1000ir_cube[k,l,*]
		    
		    if (total(check_vec EQ ignore_val) EQ 0) then begin
		    
		    	bdi1000ir_value[k,l] = cat_int_tabulated(wavelength_vec_um, 1.0 - spec_vec)
		
		    	if (keyword_set(show_plot)) then begin
			
			    wset, 0
			    plot, wvt, cube[k,l,*], $
			    	  xrange = [1000,2600], xsty=1, $
				  yrange = [0.0, 0.35], ysty=1
			    oplot, wavelength_vec, cont_cube[k,l,*]
			    wait, 0.20
			
			;stop
			
			endif    	
		
		    endif
		    
		endfor
	    endfor

    	    if (keyword_set(show_map)) then begin
    	    	window, 13, xsize = 700, ysize = 500, retain = 2
	    	tv, bytscl(bdi1000ir_value, min = 0.0, max = 0.025)
    	    endif
	    
	    pcube[*,*,0] = bdi1000ir_value
	    
	    ;print, systime(0)

   	endif



	;------------------------------------------------------------------------
	;Find IR albedo (IRA): (using R1330 to match Yves (b/c in contin, not
	;                mineralogically sensitive)
	;------------------------------------------------------------------------
	if (index_vec[1] EQ 1) then begin
		R1330 = mro_crism_lookupwv(1330,wvt)	; 395 in Zeta2
		pcube(*,*,1) = cube(*,*,R1330)
	endif

	;------------------------------------------------------------------------
	;Find olivine index (OLINDEX):
	;------------------------------------------------------------------------
	if (index_vec[2] EQ 1) then begin
		if (n_elements(R1330) EQ 0) then begin
			R1330 = mro_crism_lookupwv(1330,wvt)	; 395 in Zeta2
		endif


		R1695 = mro_crism_lookupwv(1695,wvt)	; 340 in Zeta2
		R1080 = mro_crism_lookupwv(1080,wvt)	;  in Zeta2
		R1210 = mro_crism_lookupwv(1210,wvt)	; 413 in Zeta2
		R1470 = mro_crism_lookupwv(1470,wvt)	; 374 in Zeta2

		pcube(*,*,2) = (cube(*,*,R1695) / (0.1*cube(*,*,R1080) + 0.1*cube(*,*,R1210) + 0.4*cube(*,*,R1330) + 0.4*cube(*,*,R1470) ) ) -1
	endif

	;------------------------------------------------------------------------
	;Find LCP and HCP index (LCPINDEX; HCPINDEX):
	;------------------------------------------------------------------------
	if ((index_vec[3] EQ 1) OR (index_vec[4] EQ 1)) then begin
		if (n_elements(R1080) EQ 0) then begin
			R1080 = mro_crism_lookupwv(1080,wvt)	; 395 in Zeta2
		endif
		if (n_elements(R1330) EQ 0) then begin
			R1330 = mro_crism_lookupwv(1330,wvt)	; 395 in Zeta2
		endif
		if (n_elements(R1470) EQ 0) then begin
			R1470 = mro_crism_lookupwv(1470,wvt)	; 374 in Zeta2
		endif


		R1815 = mro_crism_lookupwv(1815,wvt)	; 322 in Zeta2
		R2067 = mro_crism_lookupwv(2067,wvt)	; 283 in Zeta2

		pcube(*,*,3) = 100. * ( ( cube(*,*,R1330)-cube(*,*,R1080) ) / ( cube(*,*,R1330)+cube(*,*,R1080) )) * ( ( cube(*,*,R1330)-cube(*,*,R1815) ) / ( cube(*,*,R1330)+cube(*,*,R1815) ) )	; LCP
		pcube(*,*,4) = 100.* ( ( cube(*,*,R1470)-cube(*,*,R1080) ) / ( cube(*,*,R1470)+cube(*,*,R1080) )) * ( ( cube(*,*,R1470)-cube(*,*,R2067) ) / ( cube(*,*,R1470)+cube(*,*,R2067) ) ) ; HCP
	endif

	;------------------------------------------------------------------------
	;Find spectral variance parameter (VAR):
	; SLM modified 26 Jun 2007 to remove 1630 and 1660 nm which may have 
	; significant artifacts due to their location at the IR zone 1 zone 2 
	; filter boundary
	;------------------------------------------------------------------------
	if (index_vec[5] EQ 1) then begin
	;From Omega:
	;Fit a line from 1-2.3 microns and find variance of observed values from
	; fit values by summing in quadrature over 26 intervening wavelengths.
	; (note: ol & px will have high variance b/c of bands--other stuff should
	; have low variance.)
	;Find reflectances that we want to fit a line to:
	;For OMEGA: VARSLOPELAMS=['0.983450','2.286900']
	; IN CRISM, the shortest IR Row is 445 with a wavelength of around 1001.6 nm
	; in Zeta2 wavelength table, so we use that
	; Closest to 2.286 microns  is row 250, 2285.95 nm
	VARSLOPELAMS=[mro_crism_lookupwv(1000,wvt,/w),mro_crism_lookupwv(2287,wvt,/w)]
	VARSLOPEWL=[mro_crism_lookupwv(1000,wvt),mro_crism_lookupwv(2287,wvt)]
	VARSLOPERS=FLTARR(NX,NY,2)
	FOR k=0,N_ELEMENTS(VARSLOPELAMS)-1 DO BEGIN
	  VARSLOPERS[*,*,k]=FLOAT(cube[*,*,varslopewl(k)])
	ENDFOR
	;Get intervening wavelengths and observed reflectances at those values:
	;26 intervening wavelengths in Omega are
	;	;LAMS=['1.026400','1.055000','1.083700','1.126800','1.213000','1.328000',$
	;      '1.328000','1.371100','1.399900','1.428600','1.471700','1.471700',$
	;      '1.514800',$
	;      '1.700700','1.814300','1.814300','1.856800',$
	;      '1.927200','1.983400','2.011300','2.067000','2.122500','2.150100',$
	;      '2.163900',$
	;     '2.205100','2.232400','2.246100']
	;In the CRISM wavetables:
		wvs = [1021, 1048, 1080, 1153, 1212,1258, 1278, 1330, 1370, 1396, 1429, 1469, 1508, $
		1561, 1693, 1752, 1812, 1878, 1930, 1983, 2009, 2068, 2121, 2141, 2167, $
		2206,2233, 2253]
		LAMS = fltarr(n_elements(wvs))
		WL = intarr(n_elements(wvs))
		OBSRS=FLTARR(NX,NY,N_ELEMENTS(wvs))
		FOR  k = 0,n_elements(wvs)-1 do begin
			LAMS(k)=mro_crism_lookupwv(wvs(k),wvt,/w)
	       		WL=mro_crism_lookupwv(wvs(k),wvt)
			OBSRS[*,*,k]=FLOAT(cube[*,*,WL])
		ENDFOR
		;Fit a line and find the variance:
		VAR=FLTARR(NX,NY)
		FOR j=0,NY-1 DO BEGIN
		  FOR i=0,NX-1 DO BEGIN
		    FIT=LINFIT(float(VARSLOPELAMS),VARSLOPERS(i,j,*))
		    PREDRS=FIT[0] + FIT[1]*LAMS
		    VAR[i,j]=TOTAL((PREDRS - OBSRS[i,j,*])^2)
		  ENDFOR
		ENDFOR
		pcube(*,*,5) = var
	endif

	;------------------------------------------------------------------------
	;Find spectral slope 1 inverse (ISLOPE1): (this is here b/c use R2530 in 1 micron bd)
	;------------------------------------------------------------------------
	if (index_vec[6] EQ 1) then begin
		if (n_elements(R1815) EQ 0) then begin
			R1815 = mro_crism_lookupwv(1815,wvt)	; 322 in Zeta2
		endif

		R2530 = mro_crism_lookupwv(2530,wvt)			; 213 in Zeta2
		W2530 = mro_crism_lookupwv(2530,wvt,/w)/1000	; 2.52965 in Zeta2
		W1815 = mro_crism_lookupwv(1815,wvt,/w)/1000	; ??????? in Zeta2
		pcube(*,*,6) = ( cube(*,*,R1815) - cube(*,*,R2530) ) /  ( W2530[0] - W1815[0] )
	endif

	;------------------------------------------------------------------------
	;Find 1.435 micron band depth (BD1435): (CO2 ice)
	;------------------------------------------------------------------------
	if (index_vec[7] EQ 1) then begin
		if (n_elements(R1470) EQ 0) then begin
			R1470 = mro_crism_lookupwv(1470,wvt)	; 374 in Zeta2
		endif


		R1430 = mro_crism_lookupwv(1430,wvt)		; 380 in Zeta2
		R1370 = mro_crism_lookupwv(1370,wvt)		; 389 in Zeta2
		WL = mro_crism_lookupwv(1370,wvt,/w)
		WC = mro_crism_lookupwv(1430,wvt,/w)
		WH = mro_crism_lookupwv(1470,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer (higher wvln)  band
		b = 1-a		      ; b gets multiplied by the shorter (lower wvln) band
		pcube(*,*,7) = 1-( cube(*,*,R1430) / ( b*cube(*,*,R1370) + a*cube(*,*,R1470) ) )
	endif

	;------------------------------------------------------------------------
	;Find 1.5 micron band depth (BD1500): (H20 ice)
	; ** Note V3.1. Previous version of eq. had R1300 instead of R1330. R1300 is not
	; ** in wave table and not in OMEGA
	; SLM modified 26 Jun 2007 to reformulated parameter derived by Frank
	; Seelos and Scott Murchie with less sensitivity to instrumental noise
	;------------------------------------------------------------------------
	if (index_vec[8] EQ 1) then begin
	    ; 1.0 - (R1558 + R1505) / (R1808+ R1367)
	    R1558 = mro_crism_lookupwv(1558,wvt)
    	    R1505 = mro_crism_lookupwv(1505,wvt)
    	    R1808 = mro_crism_lookupwv(1808,wvt)
	    R1367 = mro_crism_lookupwv(1367,wvt)
    	    pcube[*,*,8] = 1.0 - (cube[*,*,R1558] + cube[*,*,R1505]) / (cube[*,*,R1808] + cube[*,*,R1367])
;		if (n_elements(R1330) EQ 0) then begin
;			R1330 = mro_crism_lookupwv(1330,wvt)	; 395 in Zeta2
;		endif
;		if (n_elements(R1695) EQ 0) then begin
;			R1695 = mro_crism_lookupwv(1695,wvt)	; 340 in Zeta2
;		endif
;
;		R1510 = mro_crism_lookupwv(1510,wvt)		; 368 in Zeta2
;		WL = mro_crism_lookupwv(1330,wvt,/w)
;		WC = mro_crism_lookupwv(1510,wvt,/w)
;		WH = mro_crism_lookupwv(1695,wvt,/w)
;		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer ;(higher wvln) band
;		b = 1-a		      				  ; b gets multiplied by the ;shorter (lower wvln) band
;		pcube(*,*,8) = 1-( cube(*,*,R1510) / ( b*cube(*,*,R1330) + ;a*cube(*,*,R1695) ) )
	endif

	;------------------------------------------------------------------------
	;Find the ratio for CO2 H2O ice mixtures (ICER1):
	;------------------------------------------------------------------------
	if (index_vec[9] EQ 1) then begin
		if (n_elements(R1430) EQ 0) then begin
			R1430 = mro_crism_lookupwv(1430,wvt)		; 380 in Zeta2
		endif
		if (n_elements(R1510) EQ 0) then begin
			R1510 = mro_crism_lookupwv(1510,wvt)		; 368 in Zeta2
		endif

	    pcube(*,*,9) = cube(*,*,R1510) / cube(*,*,R1430)
        endif

	;------------------------------------------------------------------------
	;Find the 1.7 micron band depth (BD1750):  (wvs 1750, 1660, 1815)
	; SLM modified 26 Jun 2007 to move short wavelength shoulder from 1660 to 
	; 1557 nm to avoid artifact at IR zone 1 zone 2 filter boundary
	;------------------------------------------------------------------------
	if (index_vec[10] EQ 1) then begin
		if (n_elements(R1815) EQ 0) then begin
			R1815 = mro_crism_lookupwv(1815,wvt)	; 322 in Zeta2
		endif
		R1750 = mro_crism_lookupwv(1750,wvt)		; 331 in Zeta2
		R1557 = mro_crism_lookupwv(1557,wvt)		
		WL = mro_crism_lookupwv(1557,wvt,/w)
		WC = mro_crism_lookupwv(1750,wvt,/w)
		WH = mro_crism_lookupwv(1815,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer (higher wvln)  band
		b = 1-a		      ; b gets multiplied by the shorter (lower wvln) band
		pcube(*,*,10) = 1-( cube(*,*,R1750) / ( b*cube(*,*,R1557) + a*cube(*,*,R1815) ) )
	end

	;------------------------------------------------------------------------
	;Find the 1.9 micron band depth (BD1900): (wvs 1930, 1985, 1875, 2067)
	; SLM modified 26 Jun 2007 to reformulated parameter derived by Frank
	; Seelos and Scott Murchie with less sensitivity to instrumental noise
	; and that better ignores 2-micron ice band
	;------------------------------------------------------------------------
	if (index_vec[11] EQ 1) then begin
	    ; 1.0 - (R1973 + R1927) / (R2006 + R1874)
	    R1973 = mro_crism_lookupwv(1973,wvt)
    	    R1927 = mro_crism_lookupwv(1927,wvt)
    	    R2006 = mro_crism_lookupwv(2006,wvt)
	    R1874 = mro_crism_lookupwv(1874,wvt)
        pcube[*,*,11] = 1.0 - (cube[*,*,R1973] + cube[*,*,R1927]) / (cube[*,*,R2006] + cube[*,*,R1874])
    	endif
;		if (n_elements(R2067) EQ 0) then begin
;			R2067 = mro_crism_lookupwv(2067,wvt)	; 283 in Zeta2
;		endif
;
;		R1875 = mro_crism_lookupwv(1875,wvt)		; 312 in Zeta2
;		R1930 = mro_crism_lookupwv(1930,wvt)		; 304 in Zeta2
;		R1985 = mro_crism_lookupwv(1985,wvt)		; 296 in Zeta2
;		WL = mro_crism_lookupwv(1875,wvt,/w)
;		WC = ;(mro_crism_lookupwv(1985,wvt,/w)+mro_crism_lookupwv(1930,wvt,/w))*0.5
;		WH = mro_crism_lookupwv(2067,wvt,/w)
;		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer ;(higher wvln)  band
;		b = 1-a		      ; b gets multiplied by the shorter (lower wvln) ;band
;		pcube(*,*,11) = 1 - ( ( (cube(*,*,R1930) + cube(*,*,R1985) )*0.5 ) / ( ;a*cube(*,*,R1875) + b*cube(*,*,R2067) ) )
;	endif


	;------------------------------------------------------------------------
	;Find 2 micron band depth (BDI2000):
	;   03/24/2008	(fps) 
	;   Total rewrite of the BDI2000 calculation
	;   Uses revised continuum calculation from BDI1000IR code block
    	;   Revised ~2000 nm integrated band depth calculation
	;------------------------------------------------------------------------
	if (index_vec[12] EQ 1) then begin

    	    ;print, 'BDI2000'
	    ;print, systime(0)

    	    ;Normalize and integrate ~2 um band

    	    if (n_elements(R1660) EQ 0) then begin
    	    	R1660 = mro_crism_lookupwv(1660,wvt)
	    endif
	    if (n_elements(R1690) EQ 0) then begin
    	    	R1690 = mro_crism_lookupwv(1690,wvt)
	    endif
	    if (n_elements(R1750) EQ 0) then begin
    	    	R1750 = mro_crism_lookupwv(1750,wvt)
	    endif
	    if (n_elements(R1810) EQ 0) then begin
    	    	R1810 = mro_crism_lookupwv(1810,wvt)
	    endif
	    if (n_elements(R1875) EQ 0) then begin
    	    	R1875 = mro_crism_lookupwv(1875,wvt)
	    endif
	    
	    if (n_elements(R1930) EQ 0) then begin
    	    	R1930 = mro_crism_lookupwv(1930,wvt)
	    endif
	    if (n_elements(R1975) EQ 0) then begin
    	    	R1975 = mro_crism_lookupwv(1975,wvt)
	    endif
	    if (n_elements(R1980) EQ 0) then begin
    	    	R1980 = mro_crism_lookupwv(1980,wvt)
	    endif
	    if (n_elements(R2005) EQ 0) then begin
    	    	R2005 = mro_crism_lookupwv(2005,wvt)
	    endif
	    if (n_elements(R2065) EQ 0) then begin
    	    	R2065 = mro_crism_lookupwv(2065,wvt)
	    endif
	    
	    if (n_elements(R2120) EQ 0) then begin
    	    	R2120 = mro_crism_lookupwv(2120,wvt)
	    endif
	    if (n_elements(R2140) EQ 0) then begin
    	    	R2140 = mro_crism_lookupwv(2140,wvt)
	    endif
	    if (n_elements(R2165) EQ 0) then begin
    	    	R2165 = mro_crism_lookupwv(2165,wvt)
	    endif
    	    if (n_elements(R2205) EQ 0) then begin
    	    	R2205 = mro_crism_lookupwv(2205,wvt)
	    endif
	    if (n_elements(R2230) EQ 0) then begin
    	    	R2230 = mro_crism_lookupwv(2230,wvt)
	    endif
	    
	    if (n_elements(R2250) EQ 0) then begin
    	    	R2250 = mro_crism_lookupwv(2250,wvt)
	    endif
	    if (n_elements(R2290) EQ 0) then begin
    	    	R2290 = mro_crism_lookupwv(2290,wvt)
	    endif
	    if (n_elements(R2315) EQ 0) then begin
    	    	R2315 = mro_crism_lookupwv(2315,wvt)
	    endif
	    if (n_elements(R2330) EQ 0) then begin
    	    	R2330 = mro_crism_lookupwv(2330,wvt)
	    endif
	    if (n_elements(R2350) EQ 0) then begin
    	    	R2350 = mro_crism_lookupwv(2350,wvt)
	    endif
	    if (n_elements(R2390) EQ 0) then begin
    	    	R2390 = mro_crism_lookupwv(2390,wvt)
	    endif
	    if (n_elements(R2430) EQ 0) then begin
    	    	R2430 = mro_crism_lookupwv(2430,wvt)
	    endif
	    if (n_elements(R2455) EQ 0) then begin
    	    	R2455 = mro_crism_lookupwv(2455,wvt)
	    endif

;    	    wavelength_indx = [R1660,R1690,R1750,R1810,R1875,R1930,R1975,R1980,R2005,R2065,R2120,R2140,R2165,R2205,R2230,R2250,R2290,R2315,R2330,R2350,R2390,R2430,R2455]
;    	    wavelength_indx = [R1875,R1930,R1975,R1980,R2005,R2065,R2120,R2140,R2165,R2205,R2230,R2250,R2290,R2315,R2330,R2350,R2390,R2430,R2455]
    	    wavelength_indx = [R1875,R1930,R1975,R1980,R2005,R2065,R2120,R2140,R2165,R2205,R2230,R2250,R2290,R2315,R2330,R2350,R2390]
	    wavelength_vec = wvt[wavelength_indx]

    	    wavelength_vec_um = wavelength_vec / 1000.0

    	    ;print, '  Wavelength vector:'
    	    ;print, wavelength_vec
	    ;print, '  Index vector:'
    	    ;print, wavelength_indx

    	    bdi2000_cube = cube[*,*,wavelength_indx]
    	    
	    cont_cube = make_array(nx, ny, n_elements(wavelength_vec), value = ignore_val)

	    for k = 0, n_elements(wavelength_vec) -1 do begin
    	    	cont_cube[*,*,k] = cont_slope_frame * (wavelength_vec[k] - median2_wavelength) + median2_frame
	    endfor

    	    bdi2000_normalized_cube = bdi2000_cube / cont_cube

    	    bdi2000_value = make_array(nx, ny, value = ignore_val)

	    pcube[*,*,12] = bdi2000_value
	    
	    for k = 0, nx -1 do begin
	    	for l = 0, ny -1 do begin

    	    	    spec_vec = bdi2000_normalized_cube[k,l,*]
    	    	    check_vec = bdi2000_cube[k,l,*]
		    
		    if (total(check_vec EQ ignore_val) EQ 0) then begin
		    
		    	bdi2000_value[k,l] = cat_int_tabulated(wavelength_vec_um, 1.0 - spec_vec)
		
		    endif
		    
		endfor
	    endfor

    	    if (keyword_set(show_map)) then begin
    	    	window, 14, xsize = 700, ysize = 500, retain = 2
    	    	tv, bytscl(bdi2000_value, min = 0.0, max = 0.025)
	    endif
    	    
	    pcube[*,*,12] = bdi2000_value
	    
	    ;print, systime(0)

	endif

	;------------------------------------------------------------------------
	;Find the 2.1 micron band depth (BD2100): (shifted water band) (wvs 2120, 2140, 1930, 2250)
	;------------------------------------------------------------------------
	if (index_vec[13] EQ 1) then begin
		if (n_elements(R1930) EQ 0) then begin
			R1930 = mro_crism_lookupwv(1930,wvt)
		endif

		R2120 = mro_crism_lookupwv(2120,wvt)		; 275 in Zeta2
		R2140 = mro_crism_lookupwv(2140,wvt)		; 272 in Zeta2
		R2250 = mro_crism_lookupwv(2250,wvt)		; 255 in Zeta2
		WL = mro_crism_lookupwv(1930,wvt,/w)
		WC = (mro_crism_lookupwv(2120,wvt,/w)+mro_crism_lookupwv(2140,wvt,/w))*0.5
		WH = mro_crism_lookupwv(2250,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer (higher wvln)  band
		b = 1-a		      ; b gets multiplied by the shorter (lower wvln) band
		pcube(*,*,13) = 1-( ( (cube(*,*,R2120) + cube(*,*,R2140) )*0.5 ) / ( b*cube(*,*,R1930) + a*cube(*,*,R2250) ) )
	endif

	;------------------------------------------------------------------------
	;Find the 2.21 micron AL-OH band depth (BD2210): (wvs 2210, 2140, 2250)
	;------------------------------------------------------------------------
	if (index_vec[14] EQ 1) then begin
		if (n_elements(R2140) EQ 0) then begin
			R2140 = mro_crism_lookupwv(2140,wvt)		; 272 in Zeta2
		endif
		if (n_elements(R2250) EQ 0) then begin
			R2250 = mro_crism_lookupwv(2250,wvt)		; 255 in Zeta2
		endif

		R2210 = mro_crism_lookupwv(2210,wvt)		; 262 in Zeta2
		WL = mro_crism_lookupwv(2140,wvt,/w)
		WC = mro_crism_lookupwv(2210,wvt,/w)
		WH = mro_crism_lookupwv(2250,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer (higher wvln)  band
		b = 1-a		      				  ; b gets multiplied by the shorter (lower wvln) band
		pcube(*,*,14) = 1-( cube(*,*,R2210) / ( b*cube(*,*,R2140) + a*cube(*,*,R2250) ) )
	endif

	;------------------------------------------------------------------------
	; Find the 2.3 micron Fe-OH band depth (BD2290): (really calculated at 2.293 micron)
	; (with this formulation, this is also a good CO2 ice parameter--2.293
	;  micron band depth) (wvs 2290, 2250, 2350)
	;------------------------------------------------------------------------
	if (index_vec[15] EQ 1) then begin
		if (n_elements(R2250) EQ 0) then begin
			R2250 = mro_crism_lookupwv(2250,wvt)		; 255 in Zeta2
		endif

		R2290 = mro_crism_lookupwv(2290,wvt)		; 249 in Zeta2
		R2350 = mro_crism_lookupwv(2350,wvt)		; 240 in Zeta2

		WL = mro_crism_lookupwv(2250,wvt,/w)
		WC = mro_crism_lookupwv(2290,wvt,/w)
		WH = mro_crism_lookupwv(2350,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer (higher wvln)  band
		b = 1-a		      				  ; b gets multiplied by the shorter (lower wvln) band
		; take out the following 2 statements, bug found by C.Viviano corrected
		; by F.Morgan Oct 25 2010
		;a = 0.6
		;b = 0.4
		pcube(*,*,15) = 1-( cube(*,*,R2290) / ( b*cube(*,*,R2250) + a*cube(*,*,R2350) ) )
	endif

	;------------------------------------------------------------------------
	;Find param gauging dropoff at 2.3 micron (D2300): (wvs 2290, 2330, 2320, 2120, 2170, 2210)
	; (do a simple continuum correction by dividing observed reflectances
	;  to those fit to a line from 1.8-2.53 microns)
    ; (Modified D2300 to use 2120 rather than 2140, May 2008, MFM)
	;------------------------------------------------------------------------
	if ((index_vec[16] EQ 1) OR (index_vec[17] EQ 1)) then begin ;coef used in index 17 calc
		if (n_elements(R1815) EQ 0) then begin
			R1815 = mro_crism_lookupwv(1815,wvt)	; 322 in Zeta2
		endif
		if (n_elements(R2530) EQ 0) then begin
			R2530 = mro_crism_lookupwv(2530,wvt)
		endif

		LAMB =2320   ; closest to Omega's 2320
		R2320 = mro_crism_lookupwv(LAMB,wvt)
		W2320 = mro_crism_lookupwv(LAMB,wvt,/w)
		LAMB =2170   ; closest to Omega's 2170
		R2170 = mro_crism_lookupwv(LAMB,wvt)
		W2170 = mro_crism_lookupwv(LAMB,wvt,/w)
		R2120 = mro_crism_lookupwv(2120,wvt)
		R2290 = mro_crism_lookupwv(2290,wvt)		; 249 in Zeta2
		R2350 = mro_crism_lookupwv(2350,wvt)		;  in Zeta2
		R2330 = mro_crism_lookupwv(2330,wvt)		; in Zeta2
		R2210 = mro_crism_lookupwv(2210,wvt)		;

		A= mro_crism_lookupwv(2530,wvt,/w);  = 2529
		B= mro_crism_lookupwv(1815,wvt,/w); = 1811

		slope=(cube(*,*,R2530)-cube(*,*,R1815))/(A[0]-B[0])
		C=mro_crism_lookupwv(2290,wvt,/w)
		CR2290= cube(*,*,R1815) + SLOPE * (C[0]-B[0])
		C=mro_crism_lookupwv(2320,wvt,/w)
		CR2320= cube(*,*,R1815) + SLOPE * (C[0]-B[0])
		C=mro_crism_lookupwv(2330,wvt,/w)
		CR2330= cube(*,*,R1815) + SLOPE * (C[0]-B[0])
		C=mro_crism_lookupwv(2120,wvt,/w)
		CR2120= cube(*,*,R1815) + SLOPE * (C[0]-B[0])
		C=mro_crism_lookupwv(2170,wvt,/w)
		CR2170= cube(*,*,R1815) + SLOPE * (C[0]-B[0])
		C=mro_crism_lookupwv(2210,wvt,/w)
		CR2210= cube(*,*,R1815) + SLOPE * (C[0]-B[0])
		D2300=1-( ((cube(*,*,R2290)/CR2290) + (cube(*,*,R2320)/CR2320) + (cube(*,*,R2330)/CR2330))/((cube(*,*,R2120)/CR2120) + (cube(*,*,R2170)/CR2170) + (cube(*,*,R2210)/CR2210)) )

        pcube(*,*,16)=D2300
	endif

	;------------------------------------------------------------------------
	;Find param gauging convexity at 2.29 micron due to sulfate absorptions
        ; at 1.9/2.1 and 2/4 microns (SINDEX): (wvs 2100,2400,2290)
	;------------------------------------------------------------------------
	if (index_vec[17] EQ 1) then begin
;		R2390 = mro_crism_lookupwv(2390,wvt)		; 234 in Zeta2
;		R2430 = mro_crism_lookupwv(2430,wvt)		; 228 in Zeta2
;		C=mro_crism_lookupwv(2390,wvt,/w)
;		CR2390= cube(*,*,R2390) + SLOPE * (C[0]-B[0])
;		C=mro_crism_lookupwv(2430,wvt,/w)
;		CR2430= cube(*,*,R2430) + SLOPE * (C[0]-B[0])

;		D2400= 1-( cube(*,*,R2390)/CR2390 + cube(*,*,R2430)/CR2430 ) / ( cube(*,*,R2290)/CR2290 + cube(*,*,R2320)/CR2320 )


		R2100 = mro_crism_lookupwv(2100,wvt)
		R2400 = mro_crism_lookupwv(2400,wvt)
        R2290 = mro_crism_lookupwv(2290,wvt)

        SINDEX= 1-( (cube(*,*,R2100) + cube(*,*,R2400)) /  (2*cube(*,*,R2290)) )

		pcube(*,*,17)=SINDEX

	endif

	;------------------------------------------------------------------------
	;CO2 ice discriminator (ICER2): (wvs 2530, 2600)
	;------------------------------------------------------------------------
	if (index_vec[18] EQ 1) then begin
		if (n_elements(R2530) EQ 0) then begin
			R2530 = mro_crism_lookupwv(2530,wvt)
		endif

		R2600 = mro_crism_lookupwv(2600,wvt)		; 202 in Zeta2
		pcube(*,*,18) = cube(*,*,R2530) / cube(*,*,R2600)
	endif

	;------------------------------------------------------------------------
	;Find the carbonate overtone band depth (BDCARB): (wvs 2330, 2230, 2390, 2530, 2600)
	;------------------------------------------------------------------------
	if (index_vec[19] EQ 1) then begin
		if (n_elements(R2330) EQ 0) then begin
			R2330 = mro_crism_lookupwv(2330,wvt)		; in Zeta2
		endif
		if (n_elements(R2390) EQ 0) then begin
			R2390 = mro_crism_lookupwv(2390,wvt)
		endif
		if (n_elements(R2530) EQ 0) then begin
			R2530 = mro_crism_lookupwv(2530,wvt)
		endif
		if (n_elements(R2600) EQ 0) then begin
			R2600 = mro_crism_lookupwv(2600,wvt)
		endif

		R2230 = mro_crism_lookupwv(2230,wvt)		; 258 in Zeta2
		WL1 = mro_crism_lookupwv(2230,wvt,/w)
		WC1 = (mro_crism_lookupwv(2330,wvt,/w)+mro_crism_lookupwv(2120,wvt,/w))*0.5
		WH1 = mro_crism_lookupwv(2390,wvt,/w)
		a = (WC1(0)-WL1(0))/(WH1(0)-WL1(0))   	; a gets multipled by the longer (higher wvln)  band
		b = 1-a		      						; b gets multiplied by the shorter (lower wvln) band
		WL2 = mro_crism_lookupwv(2390,wvt,/w)
		WC2 = (mro_crism_lookupwv(2530,wvt,/w)+mro_crism_lookupwv(2120,wvt,/w))*0.5
		WH2 = mro_crism_lookupwv(2600,wvt,/w)
		c = (WC2(0)-WL2(0))/(WH2(0)-WL2(0))   	; c gets multipled by the longer (higher wvln)  band
		d = 1-c		      						; d gets multiplied by the shorter (lower wvln) band
		pcube(*,*,19) = 1-( SQRT( ( cube(*,*,R2330) / ( b*cube(*,*,R2230) + a*cube(*,*,R2390) ) ) * ( cube(*,*,R2530) / ( d*cube(*,*,R2230) + c*cube(*,*,R2600) ) ) ) )
	endif

	;------------------------------------------------------------------------
	;Find the 3 micron H2O band depth (BD3000): (wvs 3000, 2530, 2210)
	;------------------------------------------------------------------------
	if (index_vec[20] EQ 1) then begin
		if (n_elements(R2210) EQ 0) then begin
			R2210 = mro_crism_lookupwv(2210,wvt)		; 262 in Zeta2
		endif
		if (n_elements(R2530) EQ 0) then begin
			R2530 = mro_crism_lookupwv(2530,wvt)
		endif

		R3000 = mro_crism_lookupwv(3000,wvt)		; 142 in Zeta2
		pcube(*,*,20) = 1-( cube(*,*,R3000) / ( cube(*,*,R2530) * ( cube(*,*,R2530) / cube(*,*,R2210) ) ) )
	endif

	;------------------------------------------------------------------------
	;Find the 3.1 micron H2O ice band depth (BD3100): (wvs 3120, 3000, 3250)
	;------------------------------------------------------------------------
	if (index_vec[21] EQ 1) then begin
		if (n_elements(R3000) EQ 0) then begin
			R3000 = mro_crism_lookupwv(3000,wvt)
		endif

		R3120 = mro_crism_lookupwv(3120,wvt)		; 123 in Zeta2
		R3250 = mro_crism_lookupwv(3250,wvt)		; 104 in Zeta2
		WL = mro_crism_lookupwv(3000,wvt,/w)
		WC = mro_crism_lookupwv(3120,wvt,/w)
		WH = mro_crism_lookupwv(3250,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer (higher wvln)  band
		b = 1-a		      ; b gets multiplied by the shorter (lower wvln) band
		pcube(*,*,21) = 1-( cube(*,*,R3120) / ( b*cube(*,*,R3000) + a*cube(*,*,R3250) ) )
	endif

	;------------------------------------------------------------------------
	;Find the 3.2 micron H2O ice band depth (BD3200): (really calculaed at 3.32 micron)
	; (wvs 3320, 3250, 3390)
	;------------------------------------------------------------------------
	if (index_vec[22] EQ 1) then begin
		if (n_elements(R3250) EQ 0) then begin
			R3250 = mro_crism_lookupwv(3250,wvt)
		endif

		R3320 = mro_crism_lookupwv(3320,wvt)		; 93 in Zeta2
		R3390 = mro_crism_lookupwv(3390,wvt)		; 82 in Zeta2
		WL = mro_crism_lookupwv(3250,wvt,/w)
		WC = mro_crism_lookupwv(3320,wvt,/w)
		WH = mro_crism_lookupwv(3390,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer (higher wvln)  band
		b = 1-a		      ; b gets multiplied by the shorter (lower wvln) band
		pcube(*,*,22) = 1-( cube(*,*,R3320) / ( b*cube(*,*,R3250) + a*cube(*,*,R3390) ) )
	endif

	;------------------------------------------------------------------------
	;Find the 3.4 micron carbonate band depth (BD3400): (wvs 3390, 3500, 3250, 3630)
	;------------------------------------------------------------------------
	if (index_vec[23] EQ 1) then begin
		if (n_elements(R3250) EQ 0) then begin
			R3250 = mro_crism_lookupwv(3250,wvt)
		endif
		if (n_elements(R3390) EQ 0) then begin
			R3390 = mro_crism_lookupwv(3390,wvt)		; 82 in Zeta2
		endif

		R3500 = mro_crism_lookupwv(3500,wvt)		; 66 in Zeta2
		R3630 = mro_crism_lookupwv(3630,wvt)		; 46 in Zeta2
		WL = mro_crism_lookupwv(3250,wvt,/w)
		WC = (mro_crism_lookupwv(3390,wvt,/w)+mro_crism_lookupwv(3500,wvt,/w))*0.5
		WH = mro_crism_lookupwv(3630,wvt,/w)
		c = (WC(0)-WL(0))/(WH(0)-WL(0))   ; c gets multipled by the longer (higher wvln)  band
		d = 1-c		      ; d gets multiplied by the shorter (lower wvln) band
		pcube(*,*,23) = 1-( ( cube(*,*,R3390) + cube(*,*,R3500) )*0.5 / ( d*cube(*,*,R3250) + c*cube(*,*,R3630) ) )
	endif

	;------------------------------------------------------------------------
	; VER3 Addition: Carbonate index (CINDEX):
	; (used as a replacement for BD3800 because CRISM won't go out to 4050nm
        ; CRISM ROW 1 (3925 nm) is very noisy in groundcal. Subsitute row 2
	; (3919.38 nm) for better result for now. (wvs 3750, 3630, 3950)
	;------------------------------------------------------------------------
	if (index_vec[24] EQ 1) then begin
		if (n_elements(R3630) EQ 0) then begin
			R3630 = mro_crism_lookupwv(3630,wvt)		; 46 in Zeta2
		endif

		a = mro_crism_lookupwv(3630,wvt,/w)/1000 ; Want in Microns
		b = mro_crism_lookupwv(3750,wvt,/w)
		c = mro_crism_lookupwv(3920,wvt,/w)
		R3750 = mro_crism_lookupwv(3750,wvt)		; 28 in Zeta2
		R3920 = mro_crism_lookupwv(3920,wvt)		; 2 in Zeta2
		pcube(*,*,24) = ((cube(*,*,R3750) + ((cube(*,*,R3750)-cube(*,*,R3630))/(B(0)-A(0)))*(C(0)-B(0)))/ cube(*,*,R3920))-1
	endif

	;************************************************************************
	;PARAMETERS FOR IR NON-ATMOSPHERICALLY-CORRECTED CALIBRATED DATA
	;************************************************************************

	;------------------------------------------------------------------------
	;Find the 1.265 micron O2 emission band (BD1270O2): (wvs 1261, 1268, 1250, 1280)
	;------------------------------------------------------------------------
	if (index_vec[25] EQ 1) then begin
		R1250 = mro_crism_lookupwv(1250,wvt)		; 407 in Zeta2
		R1261 = mro_crism_lookupwv(1261,wvt)		; 406 in Zeta2
		R1268 = mro_crism_lookupwv(1268,wvt)		; 405 in Zeta2
		R1280 = mro_crism_lookupwv(1280,wvt)		; 403 in Zeta2
		WL = mro_crism_lookupwv(1250,wvt,/w)
		WC = (mro_crism_lookupwv(1261,wvt,/w)+mro_crism_lookupwv(1268,wvt,/w))*0.5
		WH = mro_crism_lookupwv(1280,wvt,/w)
		a = 0.5
		b = 0.5
		c = (WC(0)-WL(0))/(WH(0)-WL(0))   ; c gets multipled by the longer (higher wvln)  band
		d = 1-c		      ; d gets multiplied by the shorter (lower wvln) band
		pcube(*,*,25) = 1-( ( a*cube(*,*,R1261) + b*cube(*,*,R1268) ) / ( d*cube(*,*,R1250) + c*cube(*,*,R1280) ) )
	endif

	;------------------------------------------------------------------------
	;Find the 1.4 micron H2O band depth (BD1400H2O): (really calculated at 1.3855 micron)
	; (wvs 1370, 1400, 1330, 1510)
	;------------------------------------------------------------------------
	if (index_vec[26] EQ 1) then begin
		if (n_elements(R1330) EQ 0) then begin
			R1330 = mro_crism_lookupwv(1330,wvt)	; 395 in Zeta2
		endif
		if (n_elements(R1370) EQ 0) then begin
			R1370 = mro_crism_lookupwv(1370,wvt)
		endif
		if (n_elements(R1510) EQ 0) then begin
			R1510 = mro_crism_lookupwv(1510,wvt)		; 368 in Zeta2
		endif

		R1400 = mro_crism_lookupwv(1400,wvt)		; 385 in Zeta2
		WL = mro_crism_lookupwv(1330,wvt,/w)
		WC = (mro_crism_lookupwv(1370,wvt,/w)+mro_crism_lookupwv(1400,wvt,/w))*0.5
		WH = mro_crism_lookupwv(1510,wvt,/w)
		a = 0.5
		b = 0.5
		c = (WC(0)-WL(0))/(WH(0)-WL(0))   ; c gets multipled by the longer (higher wvln)  band
		d = 1-c		      ; d gets multiplied by the shorter (lower wvln) band
		pcube(*,*,26) = 1-( ( a*cube(*,*,R1370) + b*cube(*,*,R1400) ) / ( d*cube(*,*,R1330) + c*cube(*,*,R1510) ) )
	endif

	;------------------------------------------------------------------------
	;Find the 2 micron CO2 band depth (BD2000CO2): (wvs 2010, 1815, 2170)
	;------------------------------------------------------------------------
	if (index_vec[27] EQ 1) then begin
		if (n_elements(R1815) EQ 0) then begin
			R1815 = mro_crism_lookupwv(1815,wvt)	; 322 in Zeta2
		endif
		if (n_elements(R2170) EQ 0) then begin
			LAMB = 2170   ; closest to Omega's 2170
			R2170 = mro_crism_lookupwv(LAMB,wvt)
		endif

		R2010 = mro_crism_lookupwv(2010,wvt)		; 292 in Zeta2
		WL = mro_crism_lookupwv(1815,wvt,/w)
		WC = mro_crism_lookupwv(2010,wvt,/w)
		; Correction 10/26/2010 FM WH is LAMB=2170, not 3390
		WH = mro_crism_lookupwv(LAMB,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer (higher wvln)  band
		b = 1-a		      ; b gets multiplied by the shorter (lower wvln) band
		; Correction 10/26/2010 FM don't swap/hardcode the weights
		;a = 0.4444
		;b = 0.5556
		; Correction 10/26/2010 FM a weights 2170, b weights 1815 - not other way around
		pcube(*,*,27) = 1 - ( cube(*,*,R2010) / ( b*cube(*,*,R1815) + a*cube(*,*,R2170) ) )
	endif

	;------------------------------------------------------------------------
	;Find the 2.35 micron CO band depth (BD2350): (wvs 2320, 2330, 2350, 2290, 2430)
	;------------------------------------------------------------------------
	if (index_vec[28] EQ 1) then begin
		if (n_elements(R2320) EQ 0) then begin
			LAMB = 2320   ; closest to Omega's 2320
			R2320 = mro_crism_lookupwv(LAMB,wvt)
		endif
		if (n_elements(R2330) EQ 0) then begin
			R2330 = mro_crism_lookupwv(2330,wvt)		; in Zeta2
		endif
		if (n_elements(R2350) EQ 0) then begin
			R2350 = mro_crism_lookupwv(2350,wvt)
		endif
		if (n_elements(R2290) EQ 0) then begin
			R2290 = mro_crism_lookupwv(2290,wvt)		; 249 in Zeta2
		endif
		if (n_elements(R2430) EQ 0) then begin
			R2430 = mro_crism_lookupwv(2430,wvt)		; 228 in Zeta2
		endif

		WL = mro_crism_lookupwv(2290,wvt,/w)
		WC = (mro_crism_lookupwv(2320,wvt,/w)+mro_crism_lookupwv(2330,wvt,/w)+mro_crism_lookupwv(2350,wvt,/w))*0.33333
		WH = mro_crism_lookupwv(2430,wvt,/w)
		a = 0.33333
		b = 0.33333
		c = 0.33333
		d = (WC(0)-WL(0))/(WH(0)-WL(0))   ; d gets multipled by the longer (higher wvln)  band
		ee = 1-d		      ; ee gets multiplied by the shorter (lower wvln) band
		pcube(*,*,28) = 1-( ( a*cube(*,*,R2320) + b*cube(*,*,R2330) + c*cube(*,*,R2350) ) / ( ee*cube(*,*,R2290) + d*cube(*,*,R2430) ) )
	endif

	;------------------------------------------------------------------------
	;Find the 2.6 micron H2O band depth (BD2600): VER3 Formulation (wvs 2600, 2530, 2630)
	;------------------------------------------------------------------------
	if (index_vec[29] EQ 1) then begin
		if (n_elements(R2530) EQ 0) then begin
			R2530 = mro_crism_lookupwv(2530,wvt)
		endif
		if (n_elements(R2600) EQ 0) then begin
			R2600 = mro_crism_lookupwv(2600,wvt)
		endif

		R2630 = mro_crism_lookupwv(2630,wvt)		; 198 in Zeta2
		WL = mro_crism_lookupwv(2530,wvt,/w)
		WC = mro_crism_lookupwv(2600,wvt,/w)
		WH = mro_crism_lookupwv(2630,wvt,/w)
		a = (WC(0)-WL(0))/(WH(0)-WL(0))   ; a gets multipled by the longer (higher wvln)  band
		b = 1-a		      ; b gets multiplied by the shorter (lower wvln) band
		pcube(*,*,29) = 1 - ( cube(*,*,R2600) / ( b*cube(*,*,R2530) + a*cube(*,*,R2630) ) )
	endif

	;------------------------------------------------------------------------
	;Find the IR ratio 2 (IRR2): (R2530, R2210)
	;------------------------------------------------------------------------
	if (index_vec[30] EQ 1) then begin
		if (n_elements(R2210) EQ 0) then begin
			R2210 = mro_crism_lookupwv(2210,wvt)		; 262 in Zeta2
		endif
		if (n_elements(R2530) EQ 0) then begin
			R2530 = mro_crism_lookupwv(2530,wvt)
		endif

		pcube(*,*,30) = cube(*,*,R2530) / cube(*,*,R2210)
	endif

	;------------------------------------------------------------------------
	;Find the 2.70 micron reflectance (R2700):
	; 2.7 (row 187)is in the IR gap. Need to use row 188 (2694.31 nm) wv = R2694
	;------------------------------------------------------------------------
	if (index_vec[31] EQ 1) then begin
		R2694 = mro_crism_lookupwv(2694,wvt)		; 188 in Zeta2
		pcube(*,*,31) = cube(*,*,R2694)
	endif

	;------------------------------------------------------------------------
	;Find the 2.70 micron CO2 band depth (BD2700):	(wv 2694, 2530, 2350
	;------------------------------------------------------------------------
	if (index_vec[32] EQ 1) then begin
		if (n_elements(R2694) EQ 0) then begin
			R2694 = mro_crism_lookupwv(2694,wvt)
		endif
		if (n_elements(R2530) EQ 0) then begin
			R2530 = mro_crism_lookupwv(2530,wvt)
		endif
		if (n_elements(R2350) EQ 0) then begin
			R2350 = mro_crism_lookupwv(2350,wvt)
		endif

		pcube(*,*,32) = 1-( cube(*,*,R2694) / ( cube(*,*,R2530) * ( cube(*,*,R2530) / cube(*,*,R2350) ) ) )
	endif

	;------------------------------------------------------------------------
	;Find the IR ratio 3 (IRR3): (wvs 3750, 3500)
        ; SMP: June 27, 2007 changed from R3750/R3500 due to Todd's
        ; request to avoid detection of temp variations
	;------------------------------------------------------------------------
	if (index_vec[33] EQ 1) then begin
		if (n_elements(R3500) EQ 0) then begin
			R3500 = mro_crism_lookupwv(3500,wvt)
		endif
		if (n_elements(R3390) EQ 0) then begin
			R3390 = mro_crism_lookupwv(3390,wvt)
		endif

		pcube(*,*,33) = cube(*,*,R3500) / cube(*,*,R3390)
	endif


   if (det eq 2) then pcubej(*,*,11:*)=pcube

    ;print, systime(0)
    ir_end_time = systime(1)
    
    ;print, 'IR Elapsed Time (minutes): ' + strtrim(string((ir_end_time - ir_start_time) / 60.0),2)

endif

if (det eq 2) then pcube=pcubej


;Return array of correct (requested) dimensions
;Apply ignore value mask to result

return_cube = pcube[*,*,where(index_vec_in EQ 1, num_return_bands)]

;return_cube_size = size(return_cube)
;return_nb = return_cube_size[3]

return_cube_mask = cmreplicate(cube_mask_band, num_return_bands)
;help, return_cube_mask

return_cube_mask_hist = histogram(return_cube_mask, min = ignore_val, max = ignore_val, reverse_indices = ri_return_cube_mask_hist)
if (return_cube_mask_hist[0] GT 0) then begin
    return_cube[ri_return_cube_mask_hist[ri_return_cube_mask_hist[0]:ri_return_cube_mask_hist[1]-1]] = ignore_val
endif

;stop

return, return_cube

end
