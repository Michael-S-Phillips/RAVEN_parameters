;------------------------------------------------------------------------
;  Find the 2.26 micron band depth (BD2265):  (jarosite, gibbsite)
;  CEV, May 2012
;  (CEV 3/22/13): Shifted hi_wvl to longer value for reduced correlation
;  with other BD22XX parameters.  Changed kernel widths of low/high from 3 to 5.
;------------------------------------------------------------------------
function crism_summary_bd2265,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 2210, 2265, 2340, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 3, hi_width  = 5 ) 

end
