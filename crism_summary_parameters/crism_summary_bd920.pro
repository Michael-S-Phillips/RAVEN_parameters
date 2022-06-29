function crism_summary_bd920,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 800, 920, 984, $ ;
                    low_width=5, mid_width=3, hi_width=5, hyper=hyper, ignore_val=ignore_val ) 

end
