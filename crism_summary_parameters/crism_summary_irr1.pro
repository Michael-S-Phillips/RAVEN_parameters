function crism_summary_irr1,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_ratio(cube, wvt, 800, 997, hyper=hyper, ignore_val=ignore_val) ;changed 1020 to 997 to avoid max wv gap tolerance

end
