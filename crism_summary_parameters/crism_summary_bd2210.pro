
function crism_summary_bd2210,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 2120, 2210, 2250, $
                                        hyper=hyper, ignore_val=ignore_val, $
                                        low_width = 5, $    ; BS by SLM (9/2011)
                                        mid_width = 3, $    ; BS by SLM (9/2011)
                                        hi_width  = 5 )      ; BS by SLM (9/2011)

end
