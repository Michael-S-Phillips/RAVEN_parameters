;(CEV 4/19/13): Hydroxylated Ferric Sulfate absorption at 2.238 microns (Lichtenburg et al. 2009)

function crism_summary_bd2230,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 2210, 2235, 2252, hyper=hyper, ignore_val=ignore_val, low_width = 3, mid_width = 3, hi_width  = 3 )      ; 

end
