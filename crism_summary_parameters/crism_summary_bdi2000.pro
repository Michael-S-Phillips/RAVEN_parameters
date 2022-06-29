;------------------------------------------------------------------------
;  Find 2 micron band depth (BDI2000):
;    03/24/2008  (fps) 
;    Total rewrite of the BDI2000 calculation
;    Uses revised continuum calculation from BDI1000IR code block
;    Revised ~2000 nm integrated band depth calculation
;    05/13/2013 (cev)
;    Rewrote how multispectral wavelengths are chosen - simplier.
;------------------------------------------------------------------------
function crism_summary_bdi2000,cube,wvt,hyper=hyper, ignore_val=ignore_val, bdicont=bdicont

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
        R1660_indx = mro_crism_lookupwv(1660,wvt)
        R2390_indx = mro_crism_lookupwv(2390,wvt)
        indices = [ R1660_indx, R2390_indx]
        wavelength_indx = indgen (max(indices) - min(indices) + 1) + min(indices)
    endif else begin
        wvs = [ 1660, 1690, 1750, 1810, 1875, 1930, 1975, 1980, 2005, 2065, 2120, $
                2140, 2165, 2205, 2230, 2250, 2290, 2315, 2330, 2350, 2390, 2430, 2455]
        wavelength_indx = fltarr( n_elements(wvs) )
    for  j = 0, n_elements(wvs)-1 do begin
        wavelength_indx[j] = mro_crism_lookupwv( wvs[j], wvt)
    endfor
    endelse

    wavelength_vec = wvt[wavelength_indx]
    wavelength_vec_um = wavelength_vec / 1000.0

    ; extract the related wavelengths from the cube, replacing CRISM NaN with IEEE NaN
    bdi2000_cube = crism_sumutil_to_nan( cube[*,*,wavelength_indx], ignore_val)

    cont_cube = make_array(nx, ny, n_elements(wavelength_vec), value = ignore_val)
    for k = 0, n_elements(wavelength_vec) - 1 do begin
        cont_cube[*,*,k] = bdicont.cont_slope_frame * $
                            (wavelength_vec[k] - bdicont.median2_wavelength) + $
                            bdicont.median2_frame
    endfor

    bdi2000_normalized_cube = bdi2000_cube / cont_cube
  
    bdi2000_value = make_array(nx, ny, value = ignore_val)

    ; visit each spatial pixel computing bdi2000
    for k = 0, nx -1 do begin
        for l = 0, ny -1 do begin
            ; fully expect IEEE NaN values in this spectral vector now
            spec_vec = bdi2000_normalized_cube[k,l,*]
            bdi2000_value[k,l] = cat_int_tabulated(wavelength_vec_um, 1.0 - spec_vec)
         endfor
    endfor

    bdi2000_value = crism_sumutil_from_nan ( bdi2000_value, ignore_val)
    
    return, bdi2000_value
end
