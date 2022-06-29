function crism_sumutil_ratio, cube, wvt, num_wavelength, denom_wavelength, $
            num_width = num_width, denom_width=denom_width, $
            hyper=hyper, ignore_val=ignore_val, norestore=norestore

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract single bands from the cube and replace CRISM_NAN with IEEE NAN
    num   = crism_sumutil_single(cube, wvt, num_wavelength, kernel_width=num_width, $
                            hyper=hyper, ignore_val=ignore_val, /norestore)
    denom = crism_sumutil_single(cube, wvt, denom_wavelength, kernel_width=denom_width, $
                            hyper=hyper, ignore_val=ignore_val, /norestore)
    img = num / denom

    ; By default, replace any lingering IEEE NaN values with CRISM_NAN.  If
    ;  norestore is set, leave the IEEE NaN values in the return img
    if not(keyword_set(norestore)) then begin
        img = crism_sumutil_from_nan(img, ignore_val )
    endif

    return,img

end
