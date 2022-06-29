;------------------------------------------------------------------------
;  Find the 3.4 micron carbonate band depth (BD3400): (wvs 3390, 3500, 3250, 3630)
;------------------------------------------------------------------------
function crism_summary_bd3400,cube,wvt,hyper=hyper,ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    if keyword_set(hyper) then begin

       ; hyper spectral movitivation:
        ;
        ;
        ;  split the difference between calcium and magnesium end-member doublet position
        ;    centered at 3400, extending 50 nm in each direction
        img = crism_sumutil_band_depth( cube, wvt, 3220, 3400, 3630, $         
                    hyper=hyper, ignore_val=ignore_val, mid_width = 15) 

;;;;;;;;;;;;;;;;;;;;;;;;;


     endif else begin

        ; extract individual channels, replacing CRISM_NANs with IEEE_NaNs
        R3250 = crism_sumutil_single(cube, wvt, 3250, hyper=0, /norestore); 
        R3390 = crism_sumutil_single(cube, wvt, 3390, hyper=0, /norestore); 
        R3500 = crism_sumutil_single(cube, wvt, 3500, hyper=0, /norestore); 
        R3630 = crism_sumutil_single(cube, wvt, 2250, hyper=0, /norestore); 

        ; compute the interpolated continuum values at selected wavelengths between 1815 and 2530
        WL = (mro_crism_lookupwv(3250,wvt,/w))[0]
        WC =((mro_crism_lookupwv(3390,wvt,/w))[0] + (mro_crism_lookupwv(3500,wvt,/w))[0])*0.5
        WH = (mro_crism_lookupwv(3630,wvt,/w))[0]
        a = (WC-WL)/(WH-WL)
        b = 1-a      

        ; compute multispectral version of bd3400
        img = 1-(((R3390 + R3500)*0.5 ) / (b*R3250 + a*R3630))
        
    endelse

    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img
end

