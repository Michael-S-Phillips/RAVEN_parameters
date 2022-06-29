;------------------------------------------------------------------------
;  CEV modified (~5/2012) significantly from original LCPINDEX formulation.
;  Fundamentally different calculation than the original version of the LCPINDEX, which was prone
;  to false-positives of non-HCP bearing material.  LCPINDEX was senstivite to convexity between
;  the ~1 & ~2 micron pyroxene absorptions.  LCPINDEX2 is sensitive to only 2-micron band and
;  not as suceptible to false positives and slope effects. Hyperspectral & multispectral friendly.
;------------------------------------------------------------------------
function crism_summary_lcpindex2,cube,wvt,hyper=hyper,ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

     ; extract channels from cube replacing CRISM_NAN with IEEE_NaN
    R1690 = crism_sumutil_single(cube, wvt, 1690, hyper=hyper, /norestore, kernel_width = 7)
    R1750 = crism_sumutil_single(cube, wvt, 1750, hyper=hyper, /norestore, kernel_width = 7)
    R1810 = crism_sumutil_single(cube, wvt, 1810, hyper=hyper, /norestore, kernel_width = 7) 
    R1870 = crism_sumutil_single(cube, wvt, 1870, hyper=hyper, /norestore, kernel_width = 7)
    
    R1560 = crism_sumutil_single(cube, wvt, 1560, hyper=hyper, /norestore, kernel_width = 7) 
    R2450 = crism_sumutil_single(cube, wvt, 2450, hyper=hyper, /norestore, kernel_width = 7) 
    
    ; identify nearest CRISM wavelength
    W1690 = (mro_crism_lookupwv(1690,wvt,/w))[0]
    W1750 = (mro_crism_lookupwv(1750,wvt,/w))[0]
    W1810 = (mro_crism_lookupwv(1810,wvt,/w))[0]
    W1870 = (mro_crism_lookupwv(1870,wvt,/w))[0]
        
    W1560 = (mro_crism_lookupwv(1560,wvt,/w))[0]
    W2450 = (mro_crism_lookupwv(2450,wvt,/w))[0]
            
    ; compute the corrected reflectance interpolating 
    slope = ( R2450 - R1560 ) / ( W2450 - W1560 )
    intercept = R2450 - slope * W2450

;  weighted sum of relative differences
    Rc1690 = slope * W1690 + intercept
    Rc1750 = slope * W1750 + intercept
    Rc1810 = slope * W1810 + intercept
    Rc1870 = slope * W1870 + intercept
       
    img=((1-(R1690/Rc1690))*0.2) + ((1-(R1750/Rc1750))*0.2) + ((1-(R1810/Rc1810))*0.3) + ((1-(R1870/Rc1870))*0.3)
    
        
 ;replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    
    return,img

end
