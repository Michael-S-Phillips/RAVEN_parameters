function crism_summary_bd2100,cube,wvt,hyper=hyper,ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    if keyword_set(hyper) then begin

        ; inspired by Olivier's hyperspectral implementation consisting of a
        ; 5 point median for the center wavelength ( from 2120 to 2140)
        ; and 3 point median for the high and low wavelength values
        ;     - HWT   08/12/2011
        img = crism_sumutil_band_depth(cube, wvt, 1930, 2132, 2250, hyper=hyper, ignore_val=ignore_val, low_width = 5, mid_width = 5, hi_width  = 5 )   ; (SLM 9/2011)
        

     endif else begin

        ; extract individual channels, replacing CRISM_NANs with IEEE_NaNs
        R1930 = crism_sumutil_single(cube, wvt, 1930, hyper=0, /norestore);
        R2120 = crism_sumutil_single(cube, wvt, 2120, hyper=0, /norestore); 
        R2130 = crism_sumutil_single(cube, wvt, 2130, hyper=0, /norestore); 
        R2250 = crism_sumutil_single(cube, wvt, 2250, hyper=0, /norestore); 

        ; compute the interpolated continuum values at selected wavelengths between 1815 and 2530
        WL = (mro_crism_lookupwv(1930,wvt,/w))[0]
        WC =((mro_crism_lookupwv(2120,wvt,/w))[0] + (mro_crism_lookupwv(2130,wvt,/w))[0])*0.5
        WH = (mro_crism_lookupwv(2250,wvt,/w))[0]
        a = (WC-WL)/(WH-WL)
        b = 1-a      

          img = 1-(((R2120 + R2130)*0.5 ) / (b*R1930 + a*R2250))

        ; replace the IEEE NAN values with CRISM_NAN
        img = crism_sumutil_from_nan(img, ignore_val)

    endelse

    return,img

end
