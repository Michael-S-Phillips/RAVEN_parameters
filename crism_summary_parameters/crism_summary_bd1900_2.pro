;------------------------------------------------------------------------
;  Find the 1.9 micron band depth (BD1900): (wvs 1930, 1985, 1875, 2067)
;   SLM modified 26 Jun 2007 to reformulated parameter derived by Frank
;   Seelos and Scott Murchie with less sensitivity to instrumental noise
;   and that better ignores 2-micron ice band
;
;  Modified this to use the standard wavelength interpolation for computing
;   band depth calculation (HWT + FPS)  (Aug 12, 2011)
;
;  CEV modified May 2012 to move short wavelength from 1874 to 1850 nm
;  and the long wavelength from 2006 to 2046 nm.
;  
;  CEV modified March 2013 to implement for hyperspectral and multispectral.
;------------------------------------------------------------------------
function crism_summary_bd1900_2,cube,wvt,hyper=hyper,ignore_val=ignore_val
  
if keyword_set(hyper) then begin
      
    img=crism_sumutil_band_depth ( cube, wvt, 1850, 1930, 2067, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 5, hi_width  = 5 )

endif else begin
  
    img1= crism_sumutil_band_depth(cube, wvt, 1850, 1930, 2067, ignore_val=ignore_val )
    img2= crism_sumutil_band_depth(cube, wvt, 1850, 1985, 2067, ignore_val=ignore_val )
    img=0.5*img1+0.5*img2
  
  ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan(img, ignore_val)
    
endelse

return, img

end
