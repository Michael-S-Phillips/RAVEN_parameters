;------------------------------------------------------------------------
; Find the 1 micron band depth (BDI1000VIS):
; (wait on 1 micron band center) VIS and IR
; 03/04/2008    (fps)
;   Total rewrite of BDI1000VIS calculation to take advantage of revised RPEAK1 procedure
;   Wavelength vector ~ [833,860,892,925,951,984,1023] 
;      - Selected to resolve to same detector rows in both multi- and hyper-spectral imgs
;   Integration is performed in microns to keep results in range of simple BD parameters [0,1]
;------------------------------------------------------------------------
function crism_summary_bdi1000vis,cube,wvt,hyper=hyper,rpk=rpk, $
                ignore_val=ignore_val, rpeak=rpeak_value

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    if (n_elements(rpeak_value) eq 0 ) then begin
        rpeak_value = crism_summary_rpeak1 (cube, wvt, hyper=hyper, ignore_val=ignore_val )
    endif

    ; prepare mask 
    sz = size(cube)
    nx = sz[1]  ; spatial detector
    ny = sz[2]  ; spatial along track
    nb = sz[3]  ; spectral bands

    ; select the wavelength table indices based on hyper switch
    if keyword_set(hyper) then begin
        R833_indx  = mro_crism_lookupwv( 833,wvt)
        R1023_indx = mro_crism_lookupwv(1023,wvt)
        indices=[ R833_indx, R1023_indx]
        wavelength_indx = indgen (max(indices) - min(indices) + 1) + min(indices)
        wavelength_vec = wvt[wavelength_indx]
    endif else begin
        R833_indx = mro_crism_lookupwv(833,wvt)
        R860_indx = mro_crism_lookupwv(860,wvt)
        R892_indx = mro_crism_lookupwv(892,wvt)
        R925_indx = mro_crism_lookupwv(925,wvt)
        R951_indx = mro_crism_lookupwv(951,wvt)
        R984_indx = mro_crism_lookupwv(984,wvt)
        R1023_indx = mro_crism_lookupwv(1023,wvt)
        wavelength_indx = [ R833_indx, R860_indx, R892_indx, R925_indx, $
                            R951_indx, R984_indx, R1023_indx ]
        wavelength_vec = wvt[wavelength_indx]
    endelse

    wavelength_vec_um = wavelength_vec / 1000.0
    bdi1000vis_value = make_array(nx, ny, value = ignore_val)

    ; pick off the selected wavelengths, replacing CRISM_NaN with IEEE_NaN
    bdi1000_cube = crism_sumutil_to_nan( cube[*,*,wavelength_indx], ignore_val)
    rpeak_value_cube = make_array(size = size(bdi1000_cube), value = !values.f_nan)
    ; replicate rpeak 2d image in all bands of rpeak cube
    for j = 0, n_elements(wavelength_indx) -1 do begin
        rpeak_value_cube[*,*,j] = rpeak_value
    endfor

    bdi1000_normalized_cube = bdi1000_cube / rpeak_value_cube

    for k = 0, nx -1 do begin
        for l = 0, ny -1 do begin
            ; There may be some IEEE NaNs here.  Let them trickle through
            spec_vec = bdi1000_normalized_cube[k,l,*]
            bdi1000vis_value[k,l] = cat_int_tabulated(wavelength_vec_um, 1.0 - spec_vec)

        endfor
    endfor

    ; replace the IEEE NAN values with CRISM_NAN
    bdi1000vis_value = crism_sumutil_from_nan ( bdi1000vis_value, ignore_val)

    return,bdi1000vis_value

end
