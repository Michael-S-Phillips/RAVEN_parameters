;------------------------------------------------------------------------
;  CEV modified (~5/2012) significantly from original HCPINDEX formulation.
;  Fundamentally different calculation than the original version of the HCPINDEX, which was prone
;  to false-positives of non-HCP bearing material.  HCPINDEX was senstivite to convexity between
;  the ~1 & ~2 micron pyroxene absorptions.  HCPINDEX2 is sensitive to only 2-micron band and
;  not as suceptible to false positives and slope effects. Hyperspectral & multispectral friendly.
;------------------------------------------------------------------------
function crism_summary_hcpindex2,cube,wvt,hyper=hyper,ignore_val=ignore_val

   if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

              ; extract channels from cube replacing CRISM_NAN with IEEE_NaN
    R2120 = crism_sumutil_single(cube, wvt, 2120, hyper=hyper, /norestore, kernel_width = 5) 
    R2140 = crism_sumutil_single(cube, wvt, 2140, hyper=hyper, /norestore, kernel_width = 7)
    R2230 = crism_sumutil_single(cube, wvt, 2230, hyper=hyper, /norestore, kernel_width = 7)
    R2250 = crism_sumutil_single(cube, wvt, 2250, hyper=hyper, /norestore, kernel_width = 7) 
    R2430 = crism_sumutil_single(cube, wvt, 2430, hyper=hyper, /norestore, kernel_width = 7) 
    R2460 = crism_sumutil_single(cube, wvt, 2460, hyper=hyper, /norestore, kernel_width = 7)
    
    R1690 = crism_sumutil_single(cube, wvt, 1690, hyper=hyper, /norestore, kernel_width = 7)
    
    R1810 = crism_sumutil_single(cube, wvt, 1810, hyper=hyper, /norestore, kernel_width = 7)
    R2530 = crism_sumutil_single(cube, wvt, 2530, hyper=hyper, /norestore, kernel_width = 7) 
    
    
    ; identify nearest CRISM wavelength
    W2120 = (mro_crism_lookupwv(2120,wvt,/w))[0]
    W2140 = (mro_crism_lookupwv(2140,wvt,/w))[0]
    W2230 = (mro_crism_lookupwv(2230,wvt,/w))[0]
    W2250 = (mro_crism_lookupwv(2250,wvt,/w))[0]
    W2430 = (mro_crism_lookupwv(2430,wvt,/w))[0]
    W2460 = (mro_crism_lookupwv(2460,wvt,/w))[0]
    
    W1690 = (mro_crism_lookupwv(1690,wvt,/w))[0]
        
    W1810 = (mro_crism_lookupwv(1810,wvt,/w))[0]
    W2530 = (mro_crism_lookupwv(2530,wvt,/w))[0]
        
    ; compute the corrected reflectance interpolating 
    ;slope = ( R2530 - R1810 ) / ( W2530 - W1810 )         
    slope = ( R2530 - R1690 ) / ( W2530 - W1690 )      
    intercept = R2530 - slope * W2530

;  weighted sum of relative differences
    Rc2120 = slope * W2120 + intercept
    Rc2140 = slope * W2140 + intercept
    Rc2230 = slope * W2230 + intercept
    Rc2250 = slope * W2250 + intercept
    Rc2430 = slope * W2430 + intercept
    Rc2460 = slope * W2460 + intercept

  img=((1-(R2120/Rc2120))*0.1) + ((1-(R2140/Rc2140))*0.1) + ((1-(R2230/Rc2230))*0.15) + ((1-(R2250/Rc2250))*0.3) + ((1-(R2430/Rc2430))*0.2) + ((1-(R2460/Rc2460))*0.15)

 ;replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    
    return,img

end
        