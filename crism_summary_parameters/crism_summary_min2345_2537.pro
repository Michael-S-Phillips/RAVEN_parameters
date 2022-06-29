;------------------------------------------------------------------------
;  Find the Fe/Ca-carbonate overtone band depth (MIN2342_2537): (Fe/Ca-Carbonate)
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
;5/12 (CEV) modified to new formulation that requires both absorption features to be present in order for index values > 0.
;

function crism_summary_min2345_2537, cube, wvt, hyper=hyper, ignore_val=ignore_val

     return, crism_sumutil_band_depth_min(cube, wvt, 2250, 2345, 2430, $
                                                     2430, 2537, 2602, $
                                         hyper=hyper, ignore_val=ignore_val, $
                                         low_width1 = 5, mid_width1 = 5, hi_width1 = 5, $
                                         low_width2 = 5, mid_width2 = 5, hi_width2 = 5 )

end
