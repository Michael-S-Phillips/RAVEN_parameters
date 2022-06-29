;------------------------------------------------------------------------
;  Find the 3 micron H2O band depth (BD3000): (wvs 3000, 2530, 2210)
;------------------------------------------------------------------------
function crism_summary_bd3000,cube,wvt,hyper=hyper,ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; get single bands from the cube, replacing CRISM_NAN with IEEE NAN
    R2210 = crism_sumutil_single(cube, wvt, 2210, hyper=hyper, /norestore); 262 in Zeta2
    R2530 = crism_sumutil_single(cube, wvt, 2530, hyper=hyper, /norestore); 213 in Zeta2
    R3000 = crism_sumutil_single(cube, wvt, 3000, hyper=hyper, /norestore); 142 in Zeta2

    img = 1 - (R3000/(R2530 * ( R2530 / R2210)))

    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img
end
