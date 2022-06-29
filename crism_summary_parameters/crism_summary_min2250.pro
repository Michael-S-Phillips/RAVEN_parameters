;------------------------------------------------------------------------
;  Find the 2.21 and 2.26 micron band depth (DOUB2250): (Opal)
;  Both absorptions must be present for index values > 0.
;  CEV, May 2012 
;  (CEV 3/22/13): Shifted low_wvl to longer value for reduced correlation
;  with other BD22XX parameters.
;------------------------------------------------------------------------
function crism_summary_min2250,cube,wvt,hyper=hyper,ignore_val=ignore_val

     return, crism_sumutil_band_depth_min(cube, wvt, 2165, 2210, 2350, $
                                                     2165, 2265, 2350, $
                                         hyper=hyper, ignore_val=ignore_val, $
                                         low_width1 = 5, mid_width1 = 3, hi_width1 = 5, $
                                         low_width2 = 5, mid_width2 = 3, hi_width2 = 5 )                                        
        
end
