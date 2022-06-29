;------------------------------------------------------------------------
; VER3 Addition: Carbonate index (CINDEX):
; (used as a replacement for BD3800 because CRISM won't go out to 4050nm
; CRISM ROW 1 (3925 nm) is very noisy in groundcal. Subsitute row 2
; (3919.38 nm) for better result for now. (wvs 3750, 3630, 3950)
; CEV modified May 2012 to fundamentally change calculation to find 
; convexity at 3.6 microns due to absoprtions at 3.4 and 3.9 microns.
;------------------------------------------------------------------------
function crism_summary_cindex2,cube,wvt, hyper=hyper, ignore_val=ignore_val

  return, crism_sumutil_band_depth_invert(cube, wvt, 3450, 3610, 3875, hyper=hyper, ignore_val=ignore_val, low_width = 9, mid_width = 11, hi_width  = 7 )

end
