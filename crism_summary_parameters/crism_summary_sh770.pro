;------------------------------------------------------------------------
;  Find the 0.77 micron band 'height' (SH770): (Iron oxides)
;  7/12 (CEV) Bands selected to be more sensitive to Fe-oxides than LCP
;------------------------------------------------------------------------
function crism_summary_sh770,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth_invert( cube, wvt, 716, 775, 860, hyper=hyper, ignore_val=ignore_val, low_width=3 ) ;changed 709->716 and kernel width to 3 to avoid max wvt gap tolerance
    
end
