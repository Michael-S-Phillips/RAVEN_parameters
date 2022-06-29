function crism_summary_icer2,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_ratio(cube, wvt, 2530, 2600, hyper=hyper, ignore_val=ignore_val)

end
