;------------------------------------------------------------------------
;  Find the 0.92 micron band depth (BD920): ('Fe mineralogy')
;  CEV modified May 2012 to move short wavelength closer to maximum for ferric minerals (from 800 to 748 nm)
;  CEV modified 11 Feb 2013 to move short and long wavelength closer to maximum for ferric minerals and to increase kernel width
;------------------------------------------------------------------------

function crism_summary_bd920_2,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 807, 920, 984, hyper=hyper, ignore_val=ignore_val, low_width=5, mid_width=5, hi_width=5 ) 

end
