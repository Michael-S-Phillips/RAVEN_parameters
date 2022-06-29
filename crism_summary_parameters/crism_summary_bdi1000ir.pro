;------------------------------------------------------------------------
; Find the 1 micron band depth (BDI1000IR):
; 03/19/2008    (fps) 
;   Total rewrite of the BDI1000IR calculation
;   Revised continuum calculation - used for both BDI1000IR and BDI2000
;       Find median value at multispectral wavelengths over the interval ~[2430,2600]
;       Linear fit between the median (wavelength, value) points is the continuum
;   Revised ~1000 nm integrated band depth calculation
;       Integration of normalized reflectance values over the inteveral ~[1045,1255]
; 04/02/2008    (fps)
;   Refined continuum calculation used in BDI1000IR and BDI2000 calculations
;       Find upper quartile (75th percentile) value at multispectral wavelengths over
;           the interval ~[1330,1810]
;       Linear fit between the ~2500 nm median (wavelength, value) and ~1500 nm upper
;           quartile (wavelength, value) points is the continuum
;  2011-08-02  (hwt)
;       added bdicont named parameter that is a structure containing parameters
;           needed to compute the integrated band depth continuum.  Constituent elements
;           of the structure are cont_slope_frame, median2_wavelength, median2_frame.
;           These are produced by the utility function crism_sumutil_bdicont.pro
;         
;------------------------------------------------------------------------
function crism_summary_bdi1000ir,cube,wvt,hyper=hyper, ignore_val=ignore_val, bdicont=bdicont

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; determine cube size for later use
    sz = size(cube)
    NX = sz[1] ; x = spatial (640 detector columns, unless binned)
    NY = sz[2] ; y = spatial (# of frames along track)
    NB = sz[3] ; z = spectral # wavelengths 

    ; compute the parameters needed for the model of the 1-2 um continuum
    if not(keyword_set(bdicont)) then begin
        bdicont = crism_sumutil_bdicont( cube, wvt, hyper=hyper, ignore_val=ignore_val)
    endif

    ;Normalize and integrate ~1 um band
    if keyword_set (hyper) then begin
        R1020_indx = mro_crism_lookupwv(1020,wvt)
        R1255_indx = mro_crism_lookupwv(1255,wvt)
        indices=[R1020_indx, R1255_indx]
        wavelength_indx = indgen (max(indices) - min(indices) + 1) + min(indices)
    endif else begin
        R1020_indx = mro_crism_lookupwv(1020,wvt)
        R1045_indx = mro_crism_lookupwv(1045,wvt)
        R1080_indx = mro_crism_lookupwv(1080,wvt)
        R1150_indx = mro_crism_lookupwv(1150,wvt)
        R1210_indx = mro_crism_lookupwv(1210,wvt)
        R1250_indx = mro_crism_lookupwv(1250,wvt)
        R1255_indx = mro_crism_lookupwv(1255,wvt)
        wavelength_indx = [R1045_indx,R1150_indx,R1210_indx,R1250_indx,R1255_indx]
    endelse
    wavelength_vec = wvt[wavelength_indx]
    wavelength_vec_um = wavelength_vec / 1000.0

    ; extract the related wavelengths from the cube, replacing CRISM NaN with IEEE NaN
    bdi1000ir_cube = crism_sumutil_to_nan( cube[*,*,wavelength_indx], ignore_val)

    ; create the continuum using the computed parameters
    cont_cube = make_array(nx, ny, n_elements(wavelength_vec), value = ignore_val)
    for k = 0, n_elements(wavelength_vec) - 1 do begin
        cont_cube[*,*,k] = bdicont.cont_slope_frame * $
                            (wavelength_vec[k] - bdicont.median2_wavelength) + $
                            bdicont.median2_frame
    endfor

    bdi1000ir_normalized_cube = bdi1000ir_cube / cont_cube
    bdi1000ir_value = make_array(nx, ny, value = ignore_val)

    for k = 0, nx -1 do begin
        for l = 0, ny -1 do begin
            spec_vec = bdi1000ir_normalized_cube[k,l,*]
            check_vec = bdi1000ir_cube[k,l,*]
            if (total(check_vec EQ ignore_val) EQ 0) then begin
                bdi1000ir_value[k,l] = cat_int_tabulated(wavelength_vec_um, 1.0 - spec_vec)
            endif
        endfor
    endfor

    ; replace the IEEE NAN values with CRISM_NAN
    bdi1000ir_value = crism_sumutil_from_nan ( bdi1000ir_value, ignore_val)

    return, bdi1000ir_value
end
