;
; Find the 0.53 micron band depth (BD530):
; SLM modified 26 Jun 2007 to move long wavelength shoulder out
; of VNIR filter zone boundary, from 648 to 709 nm
;
function crism_summary_bd530, cube, wvt, hyper=hyper, ignore_val=ignore_val
    return, crism_sumutil_band_depth ( cube, wvt, 440, 530, 716, $  ;change 709->716 and kernel width to 3 to avoid max wv gap tol
                hyper=hyper, ignore_val=ignore_val, hi_width=3)
end
