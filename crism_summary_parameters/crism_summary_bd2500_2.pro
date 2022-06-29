;------------------------------------------------------------------------
;  Find the 2.5 micron band depth due to Mg-Carbonate overtone (BD2500H):  (Mg-Carbonate)
;  CEV modified May 2012 to fundamentally change calculation to a band depth centered at 2505.
;  5/15/2013 CEV modified wavelengths so can be used on multispectral data as well.
;------------------------------------------------------------------------
function crism_summary_bd2500_2, cube, wvt, hyper=hyper, ignore_val=ignore_val

    return, crism_sumutil_band_depth(cube, wvt, 2364, 2480, 2570, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 5, hi_width  = 5 )

end 
