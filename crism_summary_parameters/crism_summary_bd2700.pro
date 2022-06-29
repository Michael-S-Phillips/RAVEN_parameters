;------------------------------------------------------------------------
;  Find the 2.70 micron CO2 band depth (BD2700):  (wv 2694, 2530, 2350)
;------------------------------------------------------------------------
function crism_summary_bd2700,cube,wvt, hyper=hyper, ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract individual bands from the cube replacing CRISM_NaNs with IEEE_NaNs
    R2694 = crism_sumutil_single ( cube,  wvt, 2694, hyper=hyper, ignore_val=ignore_val, /norestore )
    R2530 = crism_sumutil_single ( cube,  wvt, 2530, hyper=hyper, ignore_val=ignore_val, /norestore )
    R2350 = crism_sumutil_single ( cube,  wvt, 2350, hyper=hyper, ignore_val=ignore_val, /norestore )

    ; compute bd2700
    img = 1.0 - ( R2694 / ( R2530 * ( R2530 / R2350 )))

    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img
end 
