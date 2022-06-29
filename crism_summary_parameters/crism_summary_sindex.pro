function crism_summary_sindex,cube,wvt,hyper=hyper,ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract channels from the cube replacing CRISM_NaN with IEEE_NaN
    R2100 = crism_sumutil_single(cube, wvt, 2100, kernel_width=5, $ ; broad feature (SLM)
                        hyper=hyper, ignore_val=ignore_val, /norestore)
    R2290 = crism_sumutil_single(cube, wvt, 2290, kernel_width=7, $ ; broad feature (SLM)
                        hyper=hyper, ignore_val=ignore_val, /norestore)
    R2400 = crism_sumutil_single(cube, wvt, 2400, kernel_width=3, $ ; potentially sharp feature (SLM)
                        hyper=hyper, ignore_val=ignore_val, /norestore)

    ; compute sindex using IEEE NaN values
    img= 1-((R2100 + R2400) / (2.0*R2290))


    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img
end
