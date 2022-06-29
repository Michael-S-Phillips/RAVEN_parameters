;------------------------------------------------------------------------
;  Find the 1.3 micron band depth (BD1300):  (wvs 1080, 1300, 1600)
;  Detection of 1.25-1.3 micron absorption associated with Fe2+ substitution in plagioclase
;  CEV 1/25/2013
;------------------------------------------------------------------------
function crism_summary_bd1300,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 1080, 1320, 1750, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 15, hi_width  = 5 ) 

end
