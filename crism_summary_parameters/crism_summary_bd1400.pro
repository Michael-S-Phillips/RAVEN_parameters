;Find the 1.4 micron band depth (BD1400): hydrated minerals

function crism_summary_bd1400,cube,wvt,hyper=hyper, ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 1330, 1395, 1467, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 3, hi_width = 5 )  ; sharp feature, but not as sharp as center wavelength

end
