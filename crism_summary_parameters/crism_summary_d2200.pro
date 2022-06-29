function crism_summary_d2200, cube, wvt, hyper=hyper, ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract individual channels, replacing CRISM_NANs with IEEE_NaNs
    R1815 = crism_sumutil_single(cube, wvt, 1815, hyper=hyper, /norestore, kernel_width = 7);
       
    R2165 = crism_sumutil_single(cube, wvt, 2165, hyper=hyper, /norestore, kernel_width = 5); 
    
    R2210 = crism_sumutil_single(cube, wvt, 2210, hyper=hyper, /norestore, kernel_width = 7);     
    R2230 = crism_sumutil_single(cube, wvt, 2230, hyper=hyper, /norestore, kernel_width = 7); 
    
    R2430 = crism_sumutil_single(cube, wvt, 2430, hyper=hyper, /norestore, kernel_width=7); 2530


    ; retrieve the CRISM wavelengths nearest the requested values
    W1815 = (mro_crism_lookupwv ( 1815, wvt, /w ))[0]
    
    W2165 = (mro_crism_lookupwv ( 2165, wvt, /w ))[0]   
    
    W2210 = (mro_crism_lookupwv ( 2210, wvt, /w ))[0]
    W2230 = (mro_crism_lookupwv ( 2230, wvt, /w ))[0]

    W2430 = (mro_crism_lookupwv ( 2430, wvt, /w ))[0]

    ; compute the interpolated continuum values at selected wavelengths between 1815 and 2530
    slope = ( R2430 - R1815 ) / ( W2430 - W1815 )

    CR2165 = R1815 + slope * ( W2165 - W1815 )
    
    CR2210 = R1815 + slope * ( W2210 - W1815 )
    CR2230 = R1815 + slope * ( W2230 - W1815 )


    ; compute d2300 with IEEE NaN values in place of CRISM NaN
    img = 1. - (((R2210/CR2210) + (R2230/CR2230)) / $
               (2*(R2165/CR2165)))
               
    ; replace the IEEE NAN values with CRISM_NAN for the return
    img = crism_sumutil_from_nan ( img, ignore_val)
    
    return,img
    
end
