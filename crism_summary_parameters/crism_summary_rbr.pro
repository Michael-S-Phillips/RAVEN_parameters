function crism_summary_rbr,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_ratio(cube, wvt, 770, 440, hyper=hyper, ignore_val=ignore_val)

end
