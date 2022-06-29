function crism_summary_icer1,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_ratio(cube, wvt, 1510, 1430, hyper=hyper, ignore_val=ignore_val)

end
