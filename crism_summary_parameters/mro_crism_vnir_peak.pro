;~11/05/2010 (fps)
;   Initial devlopment - major overhaul of RPEAK1 calculation for use in derivation of linear VNIR-->IR continuua spectra
;   	- Fit 5th order polynomial to valid CRISM hyperspectral VNIR wavelengths as provided by mro_crism_wavelength_defaults()
;   	- Calculate oversampled model using fit model parameters
;   	- Wavelength corresponding to the model maximum is peak_lambda (cf. RPEAK1), peak_rho is the maximum model reflectance
;   	- If peak_lambda is found to be toward the long end of the VNIR wavlength range, otherwise suspect, or concavity validation is requested (/concavity)
;   	    - Find local maximum (1st derivative -->0.0) with negative concavity (2nd derivative < 0.0)
;05/19/2011 (fps)
;   Added default return selections if no relevant keywords selected - /lambda, /rho
;05/23/2011 (fps)
;   Added keyword /joined - if input cube is joined, need to specify this to mro_crism_wavelength_defaults() via pass through keyword
;11/02/2011 (fps)
;   Revised keyword(s) and return structure concerning the oversampled (internal) model
;   Exploring spectral concavity between RPEAK1 wavelength and the long wavelength end of the VNIR detector
;08/01/2012 (fps)
;   Reconfigured for production use - defaults to just the RPEAK1-style calculation
;   	Removed unnecessary calls to mro_crism_image_stats()
;   	Set sampling steps to production value (2001)
;12/16/2013 (fps)
;   Cleanup for parameter code library limited release
;   	Added keyword /snapshot - activate /show_plots for convineient spatial pixel
;   	Added keywords panic_wavelength and poly_order - override defaults (900 nm; 5th order)
;

function mro_crism_vnir_peak, in_cube, in_lambda, microns = microns, poly_order = poly_order, panic_wavelength = panic_wavelength, $
    	    	    	      lambda = lambda, rho = rho, $
    	    	    	      model_spectra = model_spectra, $
			      oversampled_model = oversampled_model, $
			      wavelength_struct = wavelength_struct, $
			      show_plots = show_plots, snapshot = snapshot, $
			      concavity = concavity, $
			      joined = joined

print, systime(0), ': MRO_CRISM_VNIR_PEAK - START'

ignore_val = 65535.0

if (not(keyword_set(poly_order))) then begin
    poly_order = 5
endif

if (not(keyword_set(panic_wavelength))) then begin
    panic_wavelength = 900  ;nm
endif

;in_cube_stats = mro_crism_image_stats(in_cube, /median)
;in_cube_stats = mro_crism_image_stats(in_cube, /percentile_stats)
;stop

;lambda_struct = mro_crism_wavelength_defaults(in_lambda)
if (keyword_set(joined)) then begin
    lambda_struct = mro_crism_wavelength_defaults(in_lambda, /conservative, /joined, /split_vnir, microns = microns)
endif else begin
    lambda_struct = mro_crism_wavelength_defaults(in_lambda, /conservative, /vnir, microns = microns)
endelse

;stop

min_valid_lambda = min(lambda_struct.wavelength_vector[lambda_struct.indx])
max_valid_lambda = max(lambda_struct.wavelength_vector[lambda_struct.indx])
;print, min_valid_lambda, max_valid_lambda

;oversampled_steps = 10001
;oversampled_steps = 1001
oversampled_steps = 2001

oversampled_lambda = findgen(oversampled_steps) / (oversampled_steps -1) * (max_valid_lambda - min_valid_lambda) + min_valid_lambda

in_cube_size = size(in_cube, /struct)

if (total([keyword_set(lambda), keyword_set(rho)]) EQ 0) then begin
    lambda = 1
    rho = 1
endif

if (keyword_set(lambda)) then begin
    peak_lambda = make_array(in_cube_size.dimensions[0], in_cube_size.dimensions[1], value = ignore_val)
endif else begin
    peak_lambda = -1
endelse

if (keyword_set(rho)) then begin
    peak_rho = make_array(in_cube_size.dimensions[0], in_cube_size.dimensions[1], value = ignore_val)
endif else begin
    peak_rho = -1
endelse

if (keyword_set(model_spectra)) then begin
    model_cube = make_array(in_cube_size.dimensions[0], in_cube_size.dimensions[1], in_cube_size.dimensions[2], value = ignore_val)
endif else begin
    model_cube = -1
endelse

if (keyword_set(oversampled_model)) then begin
    oversampled_model_cube = make_array(in_cube_size.dimensions[0], in_cube_size.dimensions[1], oversampled_steps, value = ignore_val)
endif else begin
    oversampled_model_cube = -1
endelse

if (keyword_set(wavelength_struct)) then begin
    wave_struct = lambda_struct
endif else begin
    wave_struct = -1
endelse

if (keyword_set(concavity)) then begin
    max_concavity = make_array(in_cube_size.dimensions[0], in_cube_size.dimensions[1], value = ignore_val)
    max_concavity_lambda = make_array(in_cube_size.dimensions[0], in_cube_size.dimensions[1], value = ignore_val)
    mean_concavity = make_array(in_cube_size.dimensions[0], in_cube_size.dimensions[1], value = ignore_val)
endif else begin
    max_concavity = -1
    max_concavity_lambda = -1
    mean_concavity = -1
endelse

;print, '!!'

for i = 0, in_cube_size.dimensions[0] -1 do begin
    ;print, i + 1, in_cube_size.dimensions[0]

    for j = 0, in_cube_size.dimensions[1] -1 do begin

    	if (total((in_cube[i,j,lambda_struct.indx]) EQ ignore_val) EQ 0) then begin

    	    model_params = poly_fit(lambda_struct.wavelength_vector[lambda_struct.indx], in_cube[i,j,lambda_struct.indx], poly_order)
    	    oversampled_model_spectrum = poly(oversampled_lambda, model_params)
    	
	    if (keyword_set(oversampled_model)) then begin
		oversampled_model_cube[i,j,*] = oversampled_model_spectrum
	    endif

	    if (keyword_set(model_spectra)) then begin
	    	model_cube[i,j,lambda_struct.indx] = poly(lambda_struct.wavelength_vector[lambda_struct.indx], model_params)
	    endif
	        	    
	    oversampled_model_max = max(oversampled_model_spectrum, select_indx)
;	    peak_lambda[i,j] = oversampled_lambda[oversampled_indx]

    	    if ((select_indx EQ (oversampled_steps -1)) OR (oversampled_lambda[select_indx] GE panic_wavelength) OR (keyword_set(concavity))) then begin
   	    	;print, i,j, select_indx
		;stop
				
		model_params_prime = poly_fit(oversampled_lambda, deriv(oversampled_lambda, oversampled_model_spectrum), poly_order -1)
		oversampled_model_prime = poly(oversampled_lambda, model_params_prime)

    	    	oversampled_model_prime_sign = sign_ternary(oversampled_model_prime)
		oversampled_model_prime_sign_delta = oversampled_model_prime_sign - shift(oversampled_model_prime_sign, -1)
		indx = where(oversampled_model_prime_sign_delta NE 0)

		model_params_prime_prime = poly_fit(oversampled_lambda, deriv(oversampled_lambda, oversampled_model_prime), poly_order -2)
    	    	oversampled_model_prime_prime = poly(oversampled_lambda, model_params_prime_prime)		

    	    	;sanpshot visuals...
		if (keyword_set(snapshot)) then begin
    	    	    if ((i EQ in_cube_size.dimensions[0]/10) AND (j EQ in_cube_size.dimensions[1]/10)) then  begin
		    	show_plots = 1
		    endif else begin
		    	show_plots = 0
		    endelse
		endif
	
		
    	    	if (keyword_set(show_plots)) then begin
    	    	
		    window, 0
    	    	    plot, lambda_struct.wavelength_vector[lambda_struct.indx], in_cube[i,j,lambda_struct.indx], psym = -1, $
		    	  xtitle = 'Wavelength', ytitle = 'Reflectance (Offset)', $
			  title = 'Data and Model (Order: ' + strtrim(string(poly_order),2) + ')', charsize = 1.5 
    	    	    oplot, oversampled_lambda, oversampled_model_spectrum - 0.02

    	    	    for k = 0, n_elements(indx) -1 do begin
			oplot, replicate(oversampled_lambda[indx[k]],2), [-10,10], linesty = 2
		    endfor

    	    	    window, 1
    	    	    plot, oversampled_lambda, oversampled_model_prime, $
		    	  xtitle = 'Wavelength', ytitle = 'Model Prime', $
			  title = 'Model Prime', charsize = 1.5 

    	    	    for k = 0, n_elements(indx) -1 do begin
			oplot, replicate(oversampled_lambda[indx[k]],2), [-10,10], linesty = 2
		    endfor

    	    	    oplot, oversampled_lambda, replicate(0.0, oversampled_steps), linesty = 1

    	    	    window, 2
    	    	    plot, oversampled_lambda, oversampled_model_prime_prime, $
		    	  xtitle = 'Wavelength', ytitle = 'Model Prime Prime', $
			  title = 'Model Prime Prime', charsize = 1.5 

    	    	    for k = 0, n_elements(indx) -1 do begin
			oplot, replicate(oversampled_lambda[indx[k]],2), [-10,10], linesty = 2
		    endfor

    	    	    oplot, oversampled_lambda, replicate(0.0, oversampled_steps), linesty = 1

    	    	    ;stop
		endif
		
		;select_indx = indx[where((oversampled_model_prime_prime[indx] LT 0) AND (indx NE (oversampled_steps -1)))]
		select_indx = indx[where(oversampled_model_prime_prime[indx] LT 0)]

    	    	;print, i,j
		;print, select_indx
		;print, oversampled_lambda[select_indx]

    	    	if (n_elements(select_indx) NE 1) then begin				
    	    	    ;print, '!!!'
		    ;print, select_indx
		    
		    select_indx = select_indx[0]
		
		    ;print, select_indx
		    ;print, oversampled_lambda[select_indx]
		endif
    	    	
    	    	select_prime_prime_max = max(oversampled_model_prime_prime[select_indx:n_elements(oversampled_model_prime_prime)-1], maxdx)
    	    	select_prime_prime_mean = mean(oversampled_model_prime_prime[select_indx:n_elements(oversampled_model_prime_prime)-1])
		
    	    	if (keyword_set(concavity)) then begin
    	    	    max_concavity[i,j] = select_prime_prime_max
		    max_concavity_lambda[i,j] = oversampled_lambda[maxdx + select_indx]
		    mean_concavity[i,j] = select_prime_prime_mean
		endif    	    	

	    endif   ;conditional branch (/concavity)

    	    if (keyword_set(lambda)) then begin
		peak_lambda[i,j] = oversampled_lambda[select_indx]
	    endif

	    if (keyword_set(rho)) then begin
		peak_rho[i,j] = oversampled_model_spectrum[select_indx]
	    endif
    
    	endif
	
    endfor ;j
    
endfor ;i

print, systime(0), ': MRO_CRISM_VNIR_PEAK - END'

return, {peak_lambda:peak_lambda, peak_rho:peak_rho, model_cube:model_cube, lambda_struct:wave_struct, oversampled_model:{oversampled_cube:oversampled_model_cube, oversampled_lambda:oversampled_lambda}, $
    	 concavity:{max_concavity:max_concavity, max_concavity_lambda:max_concavity_lambda, mean_concavity:mean_concavity}}

END
