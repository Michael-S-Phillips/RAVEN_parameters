;------------------------------------------------------------------------
;  Find the 2.20 and 2.26 micron band depth (DOUB2200H): (Hydrated silicate)
;  Old hyperspectral derivation (2.20 & 2.26 um doublet; Leah Roach)
;  (CEV) I'll point out that this formulation is sensitive to slope effects...
;  Also, it does not require the presence of both features to give a positive
;  value (unlike DOUB2200H2). It only still exists for continuity between old and new parameters.
;------------------------------------------------------------------------
function crism_summary_doub2200h,cube,wvt,hyper=hyper,ignore_val=ignore_val

        R2205 = crism_sumutil_single(cube, wvt, 2205, hyper=hyper, /norestore, kernel_width = 3) 
        R2258 = crism_sumutil_single(cube, wvt, 2258, hyper=hyper, /norestore, kernel_width = 3) 
        R2172 = crism_sumutil_single(cube, wvt, 2172, hyper=hyper, /norestore, kernel_width = 3) 
        R2311 = crism_sumutil_single(cube, wvt, 2311, hyper=hyper, /norestore, kernel_width = 3) 
        
        img=1.0-( (R2205 + R2258)/(R2172 + R2311) )
        
    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return, img
end
