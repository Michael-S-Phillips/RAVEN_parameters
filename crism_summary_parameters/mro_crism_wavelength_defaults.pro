;11/09/2010 (fps)
;   Initial development - CRISM wavelength vector wrangling...
;   One stop shopping for CRISM detector-specific reference wavelengths, conservative and liberal flagging of 'known bad' channels, and case-specific flagging of funky channels
;   Keyword overview:
;   	spec_vector - spectral data associated w/input wavelength_vector 
;   	/vnir, /ir - specify detector associated w/input wavelength vector - will attempt auto-detect if not specified
;   	ref_wavelength - specify reference wavelength - closest channel in wavelength_vector is identified
;   	/conservative_set - use conservative list of 'known bad' channels
;   	/non_finite - flag non_finite channels
;   	/crism_nan - flag CRISM NaN channels
;   	range_check - flag channels outide the specified range (applies to spec_vector)
;   	transfer_nan - store !values.f_nan in returned wavelength_vector for all flaged channels - useful for plotting
;   	/microns - input wavelenght_vector is in um (assumes nm otherwise)
;05/11/2011 (fps)
;   Updated to allow 'joined' (full range) wavelength vectors
;   	Added keyword /joined (associated with keywords /vnir, /ir)
;   Revised default bad bands specification to reflect TRDR_DS.CAT revision 2010-11-29
;05/23/2011 (fps)
;   Added keywords /split_vnir and /split_ir for use with /joined to include all channels from complementary detector in the out index list
;05/26/2011 (fps)
;   Added keyword rgb_wavelengths and related logic - similar behavior to ref_wavelength
;   	Tags rgb_wavelengths: and rgb_indx: added to return structure
;   Returned ref_wavelength (and rgb_wavelenghts) is/are sampled from input wavelength_vector - not idealized values
;07/12/2012 (fps)
;   split_vnir and split_ir logic works for joined wavelength vectors that have alread had bad bands removed
;06/05/2013 (fps)
;   Revised auto-detection of [/vnir, /ir, /joined] to look at wavelength vector distribution rather than number of elements (was only valid for hyperspectral data previously)
;07/16/2013 (fps)
;   Switched VNIR RGB CAT ENVI default wavelengths [709.680, 598.860, 533.740] for VNIR TRU browse product default wavelengths [598.860, 527.230, 442.630]
;   	VNIR RGB CAT ENVI wavelengths accessible via new keyword /cat_rgb
;   Added two validity test for the input wavelength_vector - return (-1) if it is not defined or if the first element is equal to (-1)
;

function mro_crism_wavelength_defaults, wavelength_vector, spec_vector = spec_vector, vnir = vnir, ir = ir, $
    	    	    	    	    	joined = joined, split_vnir = split_vnir, split_ir = split_ir, $
					ref_wavelength = ref_wavelength, rgb_wavelengths = rgb_wavelengths, cat_rgb = cat_rgb, $
    	    	    	    	    	conservative_set = conservative_set, non_finite = non_finite, crism_nan = crism_nan, range_check = range_check, $
					transfer_nan = transfer_nan, $
					microns = microns

ignore_val = 65535.0
nan_value = !values.f_nan

num_wavelengths = n_elements(wavelength_vector)

if (num_wavelengths EQ 0) then begin
    return, -1
endif

if (num_wavelengths[0] EQ -1) then begin
    return, -1
endif

if (keyword_set(microns)) then begin
    wavelength_vector = wavelength_vector * 1000.0
    
    if (keyword_set(ref_wavelength)) then begin
    	ref_wavelength = ref_wavelength * 1000.0
    endif
endif


if (total([keyword_set(ir), keyword_set(vnir), keyword_set(joined)]) NE 1) then begin
    if (median(wavelength_vector) LT 1000.0) then begin
    	vnir = 1
	ir = 0
	joined = 0
    endif else begin
    	wavelength_percentiles = percentiles(wavelength_vector)
	if (abs(wavelength_percentiles[1] - wavelength_percentiles[99]) LT 3000.0) then begin
	    vnir = 0
    	    ir = 1
	    joined = 0
    	endif else begin
	    vnir = 0
    	    ir = 0
	    joined = 1
	endelse
    endelse

;    if (num_wavelengths LT 275) then begin
;    	vnir = 1
;	ir = 0
;	joined = 0
;    endif else begin
;    	if (num_wavelengths LT 450) then begin
;	    vnir = 0
;    	    ir = 1
;	    joined = 0
;    	endif else begin
;	    vnir = 0
;    	    ir = 0
;	    joined = 1
;	endelse
;    endelse
endif

;print, vnir, ir, joined

select_vector = replicate(1, num_wavelengths)

if (keyword_set(non_finite)) then begin
    non_finite_indx = where(finite(wavelength_vector) NE 1, num_non_finite)

    if (num_non_finite GT 0) then begin
    	select_vector[non_finite_indx] = 0
 
    	if (keyword_set(transfer_nan)) then begin
    	    wavelength_vector[non_finite_indx] = nan_value
    	endif
    	
    endif
    
    if (keyword_set(spec_vector)) then begin
	non_finite_indx = where(finite(spec_vector) NE 1, num_non_finite)

	if (num_non_finite GT 0) then begin
    	    select_vector[non_finite_indx] = 0
    	    
	    if (keyword_set(transfer_nan)) then begin
	    	wavelength_vector[non_finite_indx] = nan_value
	    endif
	    
	endif
	
    endif
    
endif

if (keyword_set(crism_nan)) then begin
    crism_nan_indx = where(wavelength_vector EQ ignore_val, num_crism_nan)
    if (num_crism_nan GT 0) then begin
    	select_vector[crism_nan_indx] = 0
    	
	if (keyword_set(transfer_nan)) then begin
    	    wavelength_vector[num_crism_nan] = nan_value
    	endif

    endif
    
    if (keyword_set(spec_vector)) then begin
	crism_nan_indx = where(spec_vector EQ ignore_val, num_crism_nan)

	if (num_crism_nan GT 0) then begin
    	    select_vector[crism_nan_indx] = 0

    	    if (keyword_set(transfer_nan)) then begin
	    	wavelength_vector[crism_nan_indx] = nan_value
	    endif
	    
	endif
	    	
    endif

endif

if (keyword_set(range_check)) then begin
    range_indx = where((spec_vector LT range_check[0]) OR (spec_vector GT range_check[1]), num_range_indx)

    if (num_range_indx GT 0) then begin
    	select_vector[range_indx] = 0
	
    	if (keyword_set(transfer_nan)) then begin
	    wavelength_vector[range_indx] = nan_value
	endif	
    
    endif
    
endif


if (keyword_set(vnir)) then begin
    if (not(keyword_set(ref_wavelength))) then begin
    	ref_wavelength = 770.0
    endif

    if (not(keyword_set(rgb_wavelengths))) then begin
    	if (keyword_set(cat_rgb)) then begin
    	    rgb_wavelengths = [710.0, 599.0, 534.0]
	endif else begin
	    rgb_wavelengths = [600.0, 530.0, 440.0]
	endelse
    endif

;    min_lambda_delta = min(abs(wavelength_vector - ref_wavelength), ref_indx, /nan)
    if (keyword_set(conservative_set)) then begin
    	;out_indx = where((wavelength_vector LT 410) OR (wavelength_vector GT 1023) OR ((wavelength_vector GT 635) AND (wavelength_vector LT 695)), num_out_indx, complement = in_indx)    	
	out_indx = where((wavelength_vector LT 430) OR (wavelength_vector GT 1015) OR ((wavelength_vector GT 635) AND (wavelength_vector LT 706)), num_out_indx, complement = in_indx)    	
    endif else begin
	out_indx = where((wavelength_vector LT 410) OR (wavelength_vector GT 1023) OR ((wavelength_vector GT 644) AND (wavelength_vector LT 684)), num_out_indx, complement = in_indx) 
    endelse

;stop
        
endif

	
if (keyword_set(ir)) then begin
    if (not(keyword_set(ref_wavelength))) then begin
    	ref_wavelength = 1330.0   
    endif

    if (not(keyword_set(rgb_wavelengths))) then begin
    	rgb_wavelengths = [2530.0, 1507.0, 1080.0]
    endif

    if (keyword_set(conservative_set)) then begin
    	;out_indx = where((wavelength_vector LT 1044) OR (wavelength_vector GT 3850) OR ((wavelength_vector GT 2657) AND (wavelength_vector LT 2797)), num_out_indx, complement = in_indx) 
	out_indx = where((wavelength_vector LT 1044) OR (wavelength_vector GT 3900) OR ((wavelength_vector GT 2657) AND (wavelength_vector LT 2797)), num_out_indx, complement = in_indx) 
    endif else begin
        out_indx = where((wavelength_vector LT 1021) OR (wavelength_vector GT 3924) OR ((wavelength_vector GT 2693) AND (wavelength_vector LT 2702)), num_out_indx, complement = in_indx) 
    endelse
       
endif

if (keyword_set(joined)) then begin
    if (not(keyword_set(ref_wavelength))) then begin
    	ref_wavelength = 1330.0   
    endif

    if (not(keyword_set(rgb_wavelengths))) then begin
    	rgb_wavelengths = [2530.0, 1330.0, 770.0]
    endif

    ;try to identify split indx based on wavelength crossover
    split_indx = (where(wavelength_vector - shift(wavelength_vector,-1) GT 0))[0]
        
    if ((split_indx EQ num_wavelengths -1) OR (split_indx EQ -1))then begin

    	;VNIR + IR joined data - bad bands already removed...
	;find split indx by large gap near 1000nm

    	gap_wavelength = 1000.0
    	gap_lambda_delta = min(abs(wavelength_vector - gap_wavelength), gap_indx, /nan)

    	gap_search_half_width = 10
    	gap_wavelength_vector = wavelength_vector[(gap_indx - gap_search_half_width):(gap_indx + gap_search_half_width)]
    	gap_delta_min = min(gap_wavelength_vector - shift(gap_wavelength_vector, -1), gap_delta_mindx)

    	split_indx = gap_indx - gap_search_half_width + gap_delta_mindx

    endif 

    print, 'split_indx: ', split_indx
        
    vnir_wavelength_vector = wavelength_vector[0:split_indx]
    ir_wavelength_vector = wavelength_vector[split_indx+1:*]

    if (keyword_set(conservative_set)) then begin
    	;vnir_out_indx = where((vnir_wavelength_vector LT 410) OR ((vnir_wavelength_vector GT 635) AND (vnir_wavelength_vector LT 695)) OR (vnir_wavelength_vector GT 1023), num_vnir_out_indx, complement = vnir_in_indx)
	vnir_out_indx = where((vnir_wavelength_vector LT 430) OR ((vnir_wavelength_vector GT 635) AND (vnir_wavelength_vector LT 706)) OR (vnir_wavelength_vector GT 1015), num_vnir_out_indx, complement = vnir_in_indx)
	;ir_out_indx = where((ir_wavelength_vector LT 1044) OR (ir_wavelength_vector GT 3850) OR ((ir_wavelength_vector GT 2657) AND (ir_wavelength_vector LT 2797)), num_ir_out_indx, complement = ir_in_indx) 
	ir_out_indx = where((ir_wavelength_vector LT 1044) OR (ir_wavelength_vector GT 3900) OR ((ir_wavelength_vector GT 2657) AND (ir_wavelength_vector LT 2797)), num_ir_out_indx, complement = ir_in_indx) 
    endif else begin
    	vnir_out_indx = where((vnir_wavelength_vector LT 410) OR ((vnir_wavelength_vector GT 644) AND (vnir_wavelength_vector LT 684)) OR (vnir_wavelength_vector GT 1023), num_vnir_out_indx, complement = vnir_in_indx)
    	ir_out_indx = where((ir_wavelength_vector LT 1021) OR (ir_wavelength_vector GT 3924) OR ((ir_wavelength_vector GT 2693) AND (ir_wavelength_vector LT 2702)), num_ir_out_indx, complement = ir_in_indx)
    endelse


    case (1) of

	(keyword_set(split_vnir)):begin
	    out_indx = [vnir_out_indx, lindgen(n_elements(ir_wavelength_vector)) + split_indx + 1]
	    num_out_indx = num_vnir_out_indx + n_elements(ir_wavelength_vector)
	    in_indx = vnir_in_indx
	end

	(keyword_set(split_ir)):begin
	    out_indx = [lindgen(n_elements(vnir_wavelength_vector)), ir_out_indx + split_indx + 1]
	    num_out_indx = n_elements(vnir_wavelength_vector) + num_ir_out_indx
	    in_indx = ir_in_indx + split_indx + 1
	end

	else:begin
	    out_indx = [vnir_out_indx, ir_out_indx + split_indx + 1]
	    num_out_indx = num_vnir_out_indx + num_ir_out_indx
	    in_indx = [vnir_in_indx, ir_in_indx + split_indx + 1]
	endelse

    endcase

endif

min_lambda_delta = min(abs(wavelength_vector - ref_wavelength), ref_indx, /nan)
ref_wavelength = wavelength_vector[ref_indx]

rgb_indx = replicate(-1, n_elements(rgb_wavelengths))
for l = 0, n_elements(rgb_wavelengths) -1 do begin
    min_rgb_lambda_delta = min(abs(wavelength_vector - rgb_wavelengths[l]), ref_rgb_indx, /nan)
    rgb_indx[l] = ref_rgb_indx
    rgb_wavelengths[l] = wavelength_vector[ref_rgb_indx]
endfor

if (num_out_indx GT 0) then begin
    select_vector[out_indx] = 0

    if (keyword_set(transfer_nan)) then begin
    	wavelength_vector[out_indx] = nan_value
    endif
    
endif

indx = where(select_vector EQ 1, num_indx)
outdx = where(select_vector EQ 0, num_outdx)

if (keyword_set(mircons)) then begin
    wavelength_vector = wavelength_vector / 1000.0
    ref_wavelength = ref_wavelength / 1000.0
endif

return, {wavelength_vector:wavelength_vector, outdx:outdx, indx:indx, ref_wavelength:ref_wavelength, ref_indx:ref_indx, rgb_wavelengths:rgb_wavelengths, rgb_indx:rgb_indx}

END
