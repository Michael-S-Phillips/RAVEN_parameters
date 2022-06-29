function crism_summary_irr3,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_ratio(cube, wvt, 3500, 3390, $
                num_width=7, denom_width=7, $       ; noise reduction
                hyper=hyper, ignore_val=ignore_val)

end
