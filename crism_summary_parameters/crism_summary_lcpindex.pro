function crism_summary_lcpindex,cube,wvt,hyper=hyper,ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract channels from cube replacing CRISM_NAN with IEEE_NaN
    R1080 = crism_sumutil_single(cube, wvt, 1080, hyper=hyper, /norestore, kernel_width=7)
    R1330 = crism_sumutil_single(cube, wvt, 1330, hyper=hyper, /norestore, kernel_width=7) 
    R1815 = crism_sumutil_single(cube, wvt, 1815, hyper=hyper, /norestore, kernel_width=7) 

    ; compute lcp index
    img = 100. * (( R1330 - R1080 ) / ( R1330 + R1080 )) * $
                 (( R1330 - R1815 ) / ( R1330 + R1815 ))
        
            ; replace the IEEE NAN values with CRISM_NAN
        img = crism_sumutil_from_nan( img, ignore_val)
        
    return,img

end
