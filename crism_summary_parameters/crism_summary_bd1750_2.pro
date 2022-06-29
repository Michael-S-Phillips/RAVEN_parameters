;------------------------------------------------------------------------
;  Find the 1.7 micron band depth (BD1750):  (wvs 1750, 1660, 1815)
;  SLM modified 26 Jun 2007 to move short wavelength shoulder from 1660 to 
;  1557 nm to avoid artifact at IR zone 1 zone 2 filter boundary
;  CEV modified May 2012 to move short wavelength shoulder from 1553 (back) to 1690 nm
;  where the actual shoulder occurs.  This is possible as the artifact between
;  filters is greatly diminished for MTRDR products. Wide low_width is set to avoid any 
;  reminant artifacts.
;------------------------------------------------------------------------
function crism_summary_bd1750_2,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 1690, 1750, 1815, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 3, hi_width  = 5 ) 

end
