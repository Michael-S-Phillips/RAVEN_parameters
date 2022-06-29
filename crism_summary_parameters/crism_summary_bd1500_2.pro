;------------------------------------------------------------------------
; Find the combined 1.50 and 1.55 micron band depths for H2O ice
; CEV modified May 2012 to calculate an averaged band depth
; CEV modified March 2013 to make multi/hyper-friendly
;------------------------------------------------------------------------

function crism_summary_bd1500_2,cube,wvt,hyper=hyper,ignore_val=ignore_val

if keyword_set(hyper) then begin

  img = crism_sumutil_band_depth(cube, wvt, 1367, 1525, 1808, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 11, hi_width  = 5 )

endif else begin

  img1= crism_sumutil_band_depth(cube, wvt, 1367, 1505, 1808, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 7, hi_width  = 5 )
  img2= crism_sumutil_band_depth(cube, wvt, 1367, 1558, 1808, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 7, hi_width  = 5 )
  img=0.5*img1+0.5*img2
  
  ; replace the IEEE NAN values with CRISM_NAN
  img = crism_sumutil_from_nan(img, ignore_val)
  
endelse 

    return,img

end 