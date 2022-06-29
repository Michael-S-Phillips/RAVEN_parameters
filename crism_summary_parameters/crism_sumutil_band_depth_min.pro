;------------------------------------------------------------------------
;
;  This function is used to compute generic minimum band depth.
;  
;  For each of the low, middle or high wavelengths, specificy these:
;
;       low_width:  width of the kernel (spectrally) for average or median
;       low_avg:    use boxcar averaging rather than default of median
;
;  Christina Viviano
;  JHU/APL
;  05/16/2013
;
;------------------------------------------------------------------------
function crism_sumutil_band_depth_min,cube,wvt,$
            low1,mid1,hi1,$
            low2,mid2,hi2, $
            low_width1=low_width1, low_avg1=low_avg1, $
            mid_width1=mid_width1, mid_avg1=mid_avg1, $
            hi_width1=hi_width1,   hi_avg1=hi_avg1,   $
            low_width2=low_width2, low_avg2=low_avg2, $
            mid_width2=mid_width2, mid_avg2=mid_avg2, $
            hi_width2=hi_width2,   hi_avg2=hi_avg2,   $
            hyper=hyper, ignore_val=ignore_val
               
    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

img1=crism_sumutil_band_depth(cube, wvt, low1, mid1, hi1, $
     low_width=low_width1, mid_width=mid_width1, hi_width=hi_width1, $
     low_avg=low_avg1, mid_avg=mid_avg1, hi_avg=hi_avg1, $
     hyper=hyper, ignore_val=ignore_val)
     
img2=crism_sumutil_band_depth(cube, wvt, low2, mid2, hi2, $
     low_width=low_width2, mid_width=mid_width2, hi_width=hi_width2, $
     low_avg=low_avg2, mid_avg=mid_avg2, hi_avg=hi_avg2, $
     hyper=hyper, ignore_val=igmore_val)

     sz = size(cube)
        NX = sz[1]  ; spatial detector
        NY = sz[2]  ; spatial along track

        img=make_array(NX, NY, value=ignore_val)
        for i=0, NX-1 do begin
            for j=0, NY-1 do begin            
            d1=img1[i,j]
            d2=img2[i,j]
            d=[d1,d2]
            img[i,j]=min(d)
            endfor
        endfor
     
        ; replace the IEEE NAN values with ignore value.
    img = crism_sumutil_from_nan ( img, ignore_val)

    return, img


end
