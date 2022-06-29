;------------------------------------------------------------------------
;  Find the 2.1656 micron band depth (bd2165):  (Kaolinite-group)
;  CEV, 3/22/13
;------------------------------------------------------------------------
function crism_summary_bd2165,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 2120, 2165, 2230, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 3, hi_width  = 3 ) 

end
