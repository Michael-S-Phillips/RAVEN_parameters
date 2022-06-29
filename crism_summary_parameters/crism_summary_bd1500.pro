;------------------------------------------------------------------------
; Find the combined 1.50 and 1.55 micron band depths for H2O ice
;------------------------------------------------------------------------

function crism_summary_bd1500,cube,wvt,hyper=hyper, ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif
   
    ; extract individual bands with CRISM_NAN replaced with IEEE NAN
    R1367 = crism_sumutil_single(cube, wvt, 1367, hyper=hyper, ignore_val=ignore_val, /norestore)
    R1505 = crism_sumutil_single(cube, wvt, 1505, hyper=hyper, ignore_val=ignore_val, /norestore)
    R1558 = crism_sumutil_single(cube, wvt, 1558, hyper=hyper, ignore_val=ignore_val, /norestore)
    R1808 = crism_sumutil_single(cube, wvt, 1808, hyper=hyper, ignore_val=ignore_val, /norestore)

    img = 1.0 - ((R1558 + R1505) / (R1808 + R1367))

    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img
end 
