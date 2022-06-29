function crism_summary_islope1,cube,wvt,hyper=hyper,ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract individual bands with CRISM_NAN replaced with IEEE NaN
    R1815 = crism_sumutil_single( cube, wvt, 1815, hyper=hyper, $
                                ignore_val=ignore_val, /norestore )
    R2530 = crism_sumutil_single( cube, wvt, 2530, hyper=hyper, $
                                ignore_val=ignore_val, /norestore )

    W1815 = (mro_crism_lookupwv(1815,wvt,/w))[0]
    W2530 = (mro_crism_lookupwv(2530,wvt,/w))[0]

    ; want in units of reflectance / um
    img = 1000.0 * ( R1815 - R2530 )/ ( W2530 - W1815 )

    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img

end
