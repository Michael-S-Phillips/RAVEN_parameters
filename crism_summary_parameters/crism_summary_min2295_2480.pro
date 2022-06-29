;------------------------------------------------------------------------
;  Find the Mg-carbonate overtone band depth (MIN2295_2480): (Mg-Carbonate)
;  only useful for hyperspectral data
;------------------------------------------------------------------------
;
;12/30/2011 (fps)
;   It looks like there has been a long-standing error in the BDCARB parameter w/r/t the center wavelengths (WC#) used in the 
;   	calculation of the band depth weighting parameters (a,b,c,d). We can't find a way to justify the use of a wavelength (2120 nm) 
;   	that is short of both of the relevant absorption bands - the inference is that the center wavelengths (WC1, WC2) were intended 
;   	to be the average of two bands both within the feature - possibly (2330 nm and *2320*) and (2530 and *2520*). Shifting the
;   	effective center wavelengths to slightly shorter values should also make the parameter more sensitive to Mg-carbonates.
;01/03/2011 (fps)
;   Based on minimal testing the changes described above appropriately 'recenter' the parameter distribution so a featureless 
;   	spectrum has a parameter value of zero. More subtle changes to the shape of the parameter distribution have not been characterized.
;CEV modified May 2012 to new formulation that requires both absorption features to be present in order for index values > 0.
;CEV modified 4/10/2013: Only as strong as the weakest link.  Formulation fundamental change that preserves full distribution.  
;     Reports the smaller of the two band depths (that way by definition the other band exists...and is larger).
;

function crism_summary_min2295_2480,cube,wvt, hyper=hyper, ignore_val=ignore_val
     
     return, crism_sumutil_band_depth_min(cube, wvt, 2165, 2295, 2364, $
                                                     2364, 2480, 2570, $
                                         hyper=hyper, ignore_val=ignore_val, $
                                         low_width1 = 5, mid_width1 = 5, hi_width1 = 5, $
                                         low_width2 = 5, mid_width2 = 5, hi_width2 = 5 )

end
