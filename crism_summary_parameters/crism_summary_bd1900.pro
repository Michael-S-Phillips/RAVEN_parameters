;------------------------------------------------------------------------
;  Find the 1.9 micron band depth (BD1900): (wvs 1930, 1985, 1875, 2067)
;   SLM modified 26 Jun 2007 to reformulated parameter derived by Frank
;   Seelos and Scott Murchie with less sensitivity to instrumental noise
;   and that better ignores 2-micron ice band
;
;  Modified this to use the standard wavelength interpolation for computing
;   band depth calculation (HWT + FPS)  (Aug 12, 2011)
;
;------------------------------------------------------------------------
function crism_summary_bd1900,cube,wvt,hyper=hyper,ignore_val=ignore_val
    
    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

  img1= crism_sumutil_band_depth(cube, wvt, 1875, 1930, 2067, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 3, hi_width  = 5 )
  img2= crism_sumutil_band_depth(cube, wvt, 1875, 1985, 2067, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 3, hi_width  = 5 )
  img=0.5*img1+0.5*img2
  
return, img

end
