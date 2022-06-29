function crism_summary_irr2,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_ratio(cube, wvt, 2530, 2210, hyper=hyper, ignore_val=ignore_val)

end
