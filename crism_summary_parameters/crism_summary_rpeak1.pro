;
;  Compute RPEAK1 using Frank's routine
;10/18/2011 (fps)
;   Modified handling of /joined keyword to allowo calculation of RPEAK1 hyperspectral parameter from VNIR-only input 
;01/05/2011 (fps)
;   Repair multispectral branch
;

function crism_summary_rpeak1,cube,wvt,hyper=hyper,ignore_val=ignore_val, joined = joined

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ;------------------------------------------------------------------------
    ; Find reflectance peak 1 (RPEAK1):
    ; 03/04/2008    (fps) 
    ;   Total rewrite of the RPEAK1 calculation
    ;   Simple poyynomial fit of reflectances over specified wavelength range
    ;   Wavelength vector ~ [442,533,600,710,740,775,800,833,860,892] nm -
    ;       Selected to resolve to same detector rows in both multi- and hyper-spectral cubes
    ;   Bypassing ~ 660nm filter boundary as much as possible
    ;   Polynomial order reduced from 5 to 4 - compromise between bining of results
    ;        (e.g. spectral smile effect) and representation of spectral shape
    ;   Results are much improved, but RPEAK1 run time is ~4 minutes for an
    ;       FRT on crism-lnx## machine
    ; 03/31/2008  (fps)
    ;   Dramatic speed improvement using IDL poly function (improved loop data handling)
    ;------------------------------------------------------------------------
    if keyword_set(hyper) then begin

        ;
        ;  Call Frank Seelos' hyper spectral code from his catnip development tree
        ;
        print, "Calling Frank's hyper version mro_crism_vnir_peak()"
        peak_struct = mro_crism_vnir_peak(cube, wvt, joined = joined)
        rpeak_wavelength = peak_struct.peak_lambda
    	rpeak_reflectance = peak_struct.peak_rho
    	
	;stop

    endif else begin

        sz = size(cube)
        NX = sz[1]      ; x = spatial (640 detector columns, unless binned)
        NY = sz[2]      ; y = spatial (# of frames along track)
        NB = sz[3]      ; z = spectral # wavelengths 


        ; select specific wavelengths
        R442_indx = mro_crism_lookupwv(442,wvt)
        R533_indx = mro_crism_lookupwv(533,wvt)
        R600_indx = mro_crism_lookupwv(600,wvt)
        R710_indx = mro_crism_lookupwv(710,wvt)
        R740_indx = mro_crism_lookupwv(740,wvt)
        R775_indx = mro_crism_lookupwv(775,wvt)
        R800_indx = mro_crism_lookupwv(800,wvt)
        R833_indx = mro_crism_lookupwv(833,wvt)
        R860_indx = mro_crism_lookupwv(860,wvt)
        R892_indx = mro_crism_lookupwv(892,wvt)
        R925_indx = mro_crism_lookupwv(925,wvt)

        wavelength_indx = [ R442_indx, R533_indx, R600_indx, R710_indx, $
                            R740_indx, R775_indx, R800_indx, R833_indx, $
                            R860_indx, R892_indx, R925_indx ]
        wavelength_vec = wvt[wavelength_indx]

        num_model_points = 1001
        poly_degree = 4.0
        model_wavelength_vec = findgen(num_model_points) / (num_model_points - 1.0) * $
                     (max(wavelength_vec) - min(wavelength_vec)) + min(wavelength_vec)

        ;Define results band
        rpeak_wavelength = make_array(nx,ny, value = ignore_val, /float)
        rpeak_reflectance = make_array(nx,ny, value = ignore_val, /float)

        ;Isolate appropriate bands from input_cube and replace CRISM_NaN with IEEE_NaN
        rpeak_cube = crism_sumutil_to_nan( cube[*,*,wavelength_indx], ignore_val)

        ;if (show_plot EQ 1) then begin
        ;    window, 17, xsize = 600, ysize = 600, retain = 2
        ;endif

        for k = 0, nx -1 do begin
            ;print, k, systime(0)
            for l = 0, ny -1 do begin
                spec_vec = rpeak_cube[k,l,*]

                ; remember that the spectral vector might have IEEE NaNs

                ;if (show_plot EQ 1) then begin
                ;    plot, wavelength_vec, spec_vec, psym = -1, $
                ;    yrange = [0.0, max(spec_vec) + 0.10 * max(spec_vec)], ysty=1, $
                ;    xtitle = 'Wavelength (nm)', ytitle = 'Spectral value', $
                ;    title = 'Pixel x=' + strtrim(string(k),2) + ' y=' + strtrim(string(l),2)
                ;endif
                ;rpeak_params = poly_fit(wavelength_vec_finite, spec_vec_finite, poly_degree)
		rpeak_params = poly_fit(wavelength_vec, spec_vec, poly_degree)

                ;model_spec = make_array(num_model_points, value = 0.0)
                ;for m = 0, poly_degree do begin
                ;    model_spec = model_spec + rpeak_params[m] * model_wavelength_vec^(float(m))
                ;endfor
                model_spec = poly(model_wavelength_vec, rpeak_params)

                ;if (show_plot EQ 1) then begin
                ;    oplot, model_wavelength_vec, model_spec
                ;    ;wait, 0.25
                ;endif

                model_spec_max = max(model_spec, model_spec_max_indx)
                model_spec_max_wavelength = model_wavelength_vec[model_spec_max_indx]

                rpeak_wavelength[k,l] = model_spec_max_wavelength
                rpeak_reflectance[k,l] = model_spec_max
            endfor
        endfor

        ; replace the IEEE NAN values with CRISM_NAN
        ;rpeak_value = crism_sumutil_from_nan ( rpeak_value, ignore_val)
	
        rpeak_reflectance = crism_sumutil_from_nan (rpeak_reflectance, ignore_val)
	rpeak_wavelength = crism_sumutil_from_nan (rpeak_wavelength, ignore_val)
	
    endelse

    ; rpeak1 parameter has units of microns, not nanometers
    roi = where ( rpeak_wavelength ne ignore_val, count)
    if (count gt 0 ) then begin
        rpeak_wavelength(roi) = rpeak_wavelength(roi) / 1000.0
    endif

    return, {rpeak_wavelength:rpeak_wavelength, rpeak_reflectance:rpeak_reflectance}

end
