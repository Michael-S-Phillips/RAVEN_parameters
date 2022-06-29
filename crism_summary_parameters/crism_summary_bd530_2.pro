;
; Find the 0.53 micron band depth (BD530):
; SLM modified 26 Jun 2007 to move long wavelength shoulder out
; of VNIR filter zone boundary, from 648 to 709 nm
; CEV modified 11 Feb 2013 - move long wavelength shoulder shortward of gap and increased kernel widths 
; FPS modified 11 Feb 2013 - fixed typo: lo_width --> low_width 
;

function crism_summary_bd530_2, cube, wvt, hyper=hyper, ignore_val=ignore_val
    return, crism_sumutil_band_depth ( cube, wvt, 440, 530, 614, $  ;
                hyper=hyper, ignore_val=ignore_val, low_width=5, mid_width=5, hi_width=5)
end
