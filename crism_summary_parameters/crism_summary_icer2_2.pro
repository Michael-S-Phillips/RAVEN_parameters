function crism_summary_icer2_2,cube,wvt,hyper=hyper,ignore_val=ignore_val

    ;return, crism_sumutil_ratio(cube, wvt, 2530, 2600, hyper=hyper, ignore_val=ignore_val)
 
 if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif


    ; extract channels from cube replacing CRISM_NAN with IEEE_NaN
    R2600 = crism_sumutil_single(cube, wvt, 2600, hyper=hyper, /norestore, kernel_width = 5) 
    
    R2456 = crism_sumutil_single(cube, wvt, 1750, hyper=hyper, /norestore, kernel_width = 5) 
    R2530 = crism_sumutil_single(cube, wvt, 2400, hyper=hyper, /norestore, kernel_width = 5) 
    
    ; identify nearest CRISM wavelength

    W2600 = (mro_crism_lookupwv(2600,wvt,/w))[0]    
    
    W2456 = (mro_crism_lookupwv(2456,wvt,/w))[0] 
    W2530 = (mro_crism_lookupwv(2530,wvt,/w))[0]
    
    ; compute the corrected reflectance interpolating 
    slope = ( R2530 - R2456 ) / ( W2530 - W2456 )
    intercept = R2530 - slope * W2530

    ; weighted sum of relative differences
    Rc2600 = slope * W2600 + intercept

    img=(Rc2600-R2600)/(abs(Rc2600))
    
    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img
end
