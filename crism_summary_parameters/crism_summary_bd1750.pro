;------------------------------------------------------------------------
;  Find the 1.7 micron band depth (BD1750):  (wvs 1750, 1660, 1815)
;  SLM modified 26 Jun 2007 to move short wavelength shoulder from 1660 to 
;  1557 nm to avoid artifact at IR zone 1 zone 2 filter boundary
;------------------------------------------------------------------------
function crism_summary_bd1750,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 1550, 1750, 1815, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 3, hi_width  = 5 ) 
end
