;------------------------------------------------------------------------
;  Find the inverted band depth of the continuum between the 2.1 and 2.4 micron mono- and poly-hydrated sulfate features
;  CEV modified May 2012 to utilize crism_sumutil_band_depth_invert for the formulation.
;  Fundamentally different calculation than the original version of the SINDEX, which was prone
;  to false-positives of non-sulfate bearing material.
;------------------------------------------------------------------------
function crism_summary_sindex2,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth_invert(cube, wvt, 2120, 2290, 2400, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 7, hi_width  = 3)

end
