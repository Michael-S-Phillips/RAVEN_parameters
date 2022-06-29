;------------------------------------------------------------------------
; Find the 1.90 micron H2O band depth (BD1900R): (wvs 1908-1941 (in), 1862-1875, 2112-2126)
; (CEV) I'll point out that this formulation is sensitive to slope effects...and still 
; exists for continuity with old parameters.
;------------------------------------------------------------------------

function crism_summary_bd1900r,cube,wvt,hyper=hyper,ignore_val=ignore_val

    R1908 = crism_sumutil_single(cube, wvt, 1908, hyper=hyper, /norestore, kernel_width = 1) 
    R1914 = crism_sumutil_single(cube, wvt, 1914, hyper=hyper, /norestore, kernel_width = 1) 
    R1921 = crism_sumutil_single(cube, wvt, 1921, hyper=hyper, /norestore, kernel_width = 1) 
    R1928 = crism_sumutil_single(cube, wvt, 1928, hyper=hyper, /norestore, kernel_width = 1) 
    R1934 = crism_sumutil_single(cube, wvt, 1934, hyper=hyper, /norestore, kernel_width = 1) 
    R1941 = crism_sumutil_single(cube, wvt, 1941, hyper=hyper, /norestore, kernel_width = 1) 
    
    R1862 = crism_sumutil_single(cube, wvt, 1862, hyper=hyper, /norestore, kernel_width = 1) 
    R1869 = crism_sumutil_single(cube, wvt, 1869, hyper=hyper, /norestore, kernel_width = 1) 
    R1875 = crism_sumutil_single(cube, wvt, 1875, hyper=hyper, /norestore, kernel_width = 1) 
    R2112 = crism_sumutil_single(cube, wvt, 2112, hyper=hyper, /norestore, kernel_width = 1) 
    R2120 = crism_sumutil_single(cube, wvt, 2120, hyper=hyper, /norestore, kernel_width = 1) 
    R2126 = crism_sumutil_single(cube, wvt, 2126, hyper=hyper, /norestore, kernel_width = 1) 
    
    img= 1.0-( (R1908+R1914+R1921+R1928+R1934+R1941) / $
            (R1862+R1869+R1875+R2112+R2120+R2126) )
  
  ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan(img, ignore_val)

return, img

end
