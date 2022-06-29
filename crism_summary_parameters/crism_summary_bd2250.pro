;------------------------------------------------------------------------
;  Find the 2.25 micron (broader) band depth (bd2250):  (opal)
;  CEV, July 2012
;  (CEV 3/22/13: Shifted hi_wvl to longer wvl to reduce correlation
;  with other BD22XX parameters. 
;------------------------------------------------------------------------
function crism_summary_bd2250,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 2120, 2245, 2340, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 7, hi_width  = 3 ) 

end
