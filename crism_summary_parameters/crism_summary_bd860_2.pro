;------------------------------------------------------------------------
;  Find the 0.86 micron band depth (BD860): ('hematite band')
;  SLM modified 26 Jun 2007 to move long wavelength shoulder 
;  from 920 to 984 nm to sample more spectral curvature
;  CEV modified May 2012 to move short wavelength closer to maximum for ferric minerals (from 800 to 748 nm)
;  CEV modified 11 Feb 2013 to move short wavelength slightly to increase kernel size
;  FPS modified 11 Feb 2013 - fixed typo: lo_width --> low_width 
;------------------------------------------------------------------------
function crism_summary_bd860_2,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 755, 860, 977, hyper=hyper, ignore_val=ignore_val, low_width=5, mid_width=5, hi_width=5) 
    
end
