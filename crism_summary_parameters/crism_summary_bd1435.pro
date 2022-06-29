function crism_summary_bd1435,cube,wvt,hyper=hyper, ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 1370, 1434, 1470, hyper=hyper, ignore_val=ignore_val, low_width = 3, mid_width = 1, hi_width = 3 )  ; sharp feature, but not as sharp as center wavelength

end
