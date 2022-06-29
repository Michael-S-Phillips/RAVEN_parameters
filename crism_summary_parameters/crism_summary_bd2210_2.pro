;(CEV 3/22/13): Shifted short wvl from 2120 to 2165 to differentiate from
;kaolinite.

function crism_summary_bd2210_2,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 2165, 2210, 2290, $
                                        hyper=hyper, ignore_val=ignore_val, $
                                        low_width = 5, $    ; BS by SLM (9/2011)
                                        mid_width = 5, $    ; BS by SLM (9/2011)
                                        hi_width  = 5 )      ; BS by SLM (9/2011)

end
