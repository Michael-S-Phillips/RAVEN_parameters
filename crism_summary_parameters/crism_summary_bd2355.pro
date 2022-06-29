;------------------------------------------------------------------------
;  Find the 2.355 micron band depth (BD2355): (Chlorite)
;  CEV, Jan 2013 
;------------------------------------------------------------------------
function crism_summary_bd2355,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 2300, 2355, 2450, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 5, hi_width  = 5 )

end
