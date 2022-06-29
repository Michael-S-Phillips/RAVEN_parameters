;------------------------------------------------------------------------
; Find the 0.60 micron shoulder height (SH600):
; SLM modified 26 Jun 2007 to move long wavelength shoulder out
; of VNIR filter zone boundary, from 680 to 709 
; 03/04/2008 (fps)
;  Changed shoulder edge wavelengths reference to [533,710] from [530,709] 
;  to resolve the same detector rows in multi- and hyper-spectral data 
;  Reformulated as an 'inverted' band depth 
; CEV modified May 2012 to simplify code by using crism_sumutil_band_depth_invert
;------------------------------------------------------------------------
function crism_summary_sh600_2,cube,wvt,hyper=hyper,ignore_val=ignore_val
    
    return, crism_sumutil_band_depth_invert(cube, wvt, 533, 600, 716, hyper=hyper, ignore_val=ignore_val, hi_width  = 3 ); 709->716 kernel width 3 to avoid max wvt gap tol
 
end
