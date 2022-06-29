;------------------------------------------------------------------------
; Find the 0.60 micron shoulder height (SH600):
; SLM modified 26 Jun 2007 to move long wavelength shoulder out
; of VNIR filter zone boundary, from 680 to 709 
; 03/04/2008 (fps)
;  Changed shoulder edge wavelengths reference to [533,710] from [530,709] 
;  to resolve the same detector rows in multi- and hyper-spectral data 
;  Reformulated as an 'inverted' band depth
; 05/2012 (cev)
;   Re-wrote to utilize crism_sumutil_band_depth_invert
;------------------------------------------------------------------------
function crism_summary_sh600,cube,wvt,hyper=hyper,ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    img=crism_sumutil_band_depth_invert(cube, wvt, 533, 600, 710, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 5, hi_width  = 5 )
 
    return,img
end
