function crism_summary_d2300, cube, wvt, hyper=hyper, ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract individual channels, replacing CRISM_NANs with IEEE_NaNs
    R1815 = crism_sumutil_single(cube, wvt, 1815, hyper=hyper, /norestore);
    R2120 = crism_sumutil_single(cube, wvt, 2120, hyper=hyper, /norestore); 
    R2170 = crism_sumutil_single(cube, wvt, 2170, hyper=hyper, /norestore); 
    R2210 = crism_sumutil_single(cube, wvt, 2210, hyper=hyper, /norestore); 
    R2290 = crism_sumutil_single(cube, wvt, 2290, hyper=hyper, /norestore, kernel_width=3); 
    R2320 = crism_sumutil_single(cube, wvt, 2320, hyper=hyper, /norestore, kernel_width=3); 
    R2330 = crism_sumutil_single(cube, wvt, 2330, hyper=hyper, /norestore, kernel_width=3); 
    R2530 = crism_sumutil_single(cube, wvt, 2530, hyper=hyper, /norestore); 2530


    ; retrieve the CRISM wavelengths nearest the requested values
    W1815 = (mro_crism_lookupwv ( 1815, wvt, /w ))[0]
    W2120 = (mro_crism_lookupwv ( 2120, wvt, /w ))[0]
    W2170 = (mro_crism_lookupwv ( 2170, wvt, /w ))[0]
    W2210 = (mro_crism_lookupwv ( 2210, wvt, /w ))[0]
    W2290 = (mro_crism_lookupwv ( 2290, wvt, /w ))[0]
    W2320 = (mro_crism_lookupwv ( 2320, wvt, /w ))[0]
    W2330 = (mro_crism_lookupwv ( 2330, wvt, /w ))[0]
    W2530 = (mro_crism_lookupwv ( 2530, wvt, /w ))[0]

    ; compute the interpolated continuum values at selected wavelengths between 1815 and 2530
    slope = ( R2530 - R1815 ) / ( W2530 - W1815 )
    CR2120 = R1815 + slope * ( W2120 - W1815 )
    CR2170 = R1815 + slope * ( W2170 - W1815 )
    CR2210 = R1815 + slope * ( W2210 - W1815 )
    CR2290 = R1815 + slope * ( W2290 - W1815 )
    CR2320 = R1815 + slope * ( W2320 - W1815 )
    CR2330 = R1815 + slope * ( W2330 - W1815 )


    ; compute d2300 with IEEE NaN values in place of CRISM NaN
    img = 1 - (((R2290/CR2290) + (R2320/CR2320) + (R2330/CR2330)) / $
               ((R2120/CR2120) + (R2170/CR2170) + (R2210/CR2210))) 
 
    ; replace the IEEE NAN values with CRISM_NAN for the return
    img = crism_sumutil_from_nan ( img, ignore_val)
    
    return,img
    
end
