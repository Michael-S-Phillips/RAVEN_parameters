;------------------------------------------------------------------------
;  Find the 2.190 micron band depth (BD2190): (Beidellite, Allophane)
;  CEV, May 2012 
;------------------------------------------------------------------------
function crism_summary_bd2190,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth ( cube, wvt, 2120, 2185, 2250, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 3, hi_width  = 3 )

end
