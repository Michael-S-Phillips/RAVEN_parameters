;------------------------------------------------------------------------
;  Find the 2.16 and 2.21 micron band depth (DOUB2200): (Kaolinite group)
;  Both absorptions must be present for index values > 0.
;  CEV, May 2012 
;------------------------------------------------------------------------
function crism_summary_min2200,cube,wvt,hyper=hyper,ignore_val=ignore_val
   
     return, crism_sumutil_band_depth_min(cube, wvt, 2120, 2165, 2350, $
                                                     2120, 2210, 2350, $
                                         hyper=hyper, ignore_val=ignore_val, $
                                         low_width1 = 5, mid_width1 = 3, hi_width1 = 5, $
                                         low_width2 = 5, mid_width2 = 3, hi_width2 = 5 )
end
