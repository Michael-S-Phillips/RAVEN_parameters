function crism_summary_bd2500h, cube, wvt, hyper=hyper, ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract individual bands with CRISM_NAN replaced with IEEE NAN
    R2380 = crism_sumutil_single(cube, wvt, 2380, hyper=hyper, ignore_val=ignore_val, /norestore)
    R2500 = crism_sumutil_single(cube, wvt, 2500, hyper=hyper, ignore_val=ignore_val, /norestore)
    R2510 = crism_sumutil_single(cube, wvt, 2510, hyper=hyper, ignore_val=ignore_val, /norestore)
    R2540 = crism_sumutil_single(cube, wvt, 2540, hyper=hyper, ignore_val=ignore_val, /norestore)

    img = 1.0 - ((R2500 + R2510) / (R2540 + R2380))

    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img
end 
