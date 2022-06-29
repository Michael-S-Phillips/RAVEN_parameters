;------------------------------------------------------------------------
;  Find the carbonate overtone band depth (BDCARB): (wvs 2330, 2230, 2390, 2530, 2600)
;------------------------------------------------------------------------
;
;12/30/2011 (fps)
;   It looks like there has been a long-standing error in the BDCARB parameter w/r/t the center wavelengths (WC#) used in the 
;   	calculation of the band depth weighting parameters (a,b,c,d). We can't find a way to justify the use of a wavelength (2120 nm) 
;   	that is short of both of the relevant absorption bands - the inference is that the center wavelengths (WC1, WC2) were intended 
;   	to be the average of two bands both within the feature - possibly (2330 nm and *2320*) and (2530 and *2520*). Shifting the
;   	effective center wavelengths to slightly shorter values should also make the parameter more sensitive to Mg-carbonates.
;01/03/2011 (fps)
;   Based on minimal testing the changes described above appropriately 'recenter' the parameter distribution so a featureless 
;   	spectrum has a parameter value of zero. More subtle changes to the shape of the parameter distribution have not been characterized.
;

function crism_summary_bdcarb,cube,wvt, hyper=hyper, ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract channels, replacing CRISM_NAN with IEEE NAN
    R2230 = crism_sumutil_single( cube, wvt, 2230, hyper=hyper, ignore_val=ignore_val, /norestore )
    R2320 = crism_sumutil_single( cube, wvt, 2320, hyper=hyper, ignore_val=ignore_val, /norestore )
    R2330 = crism_sumutil_single( cube, wvt, 2330, hyper=hyper, ignore_val=ignore_val, /norestore )
    R2390 = crism_sumutil_single( cube, wvt, 2390, hyper=hyper, ignore_val=ignore_val, /norestore )
    R2520 = crism_sumutil_single( cube, wvt, 2520, hyper=hyper, ignore_val=ignore_val, /norestore )
    R2530 = crism_sumutil_single( cube, wvt, 2530, hyper=hyper, ignore_val=ignore_val, /norestore )
    R2600 = crism_sumutil_single( cube, wvt, 2600, hyper=hyper, ignore_val=ignore_val, /norestore )

    ; identify nearest CRISM wavelengths
    WL1 =  (mro_crism_lookupwv(2230,wvt,/w))[0]
    WC1 = ((mro_crism_lookupwv(2330,wvt,/w))[0] + $
           (mro_crism_lookupwv(2320,wvt,/w))[0])*0.5
    WH1 =  (mro_crism_lookupwv(2390,wvt,/w))[0]
    a =  ( WC1 - WL1 ) / ( WH1 - WL1 )  ; a gets multipled by the longer (higher wvln)  band
    b = 1.0-a                           ; b gets multiplied by the shorter (lower wvln) band

    WL2 =  (mro_crism_lookupwv(2390,wvt,/w))[0]
    WC2 = ((mro_crism_lookupwv(2530,wvt,/w))[0] + $
           (mro_crism_lookupwv(2520,wvt,/w))[0])*0.5
    WH2 =  (mro_crism_lookupwv(2600,wvt,/w))[0]
    c = ( WC2 - WL2 ) / ( WH2 - WL2 )   ; c gets multipled by the longer (higher wvln)  band
    d = 1.0-c                           ; d gets multiplied by the shorter (lower wvln) band

    ; compute bdcarb
    img = 1.0 - ( SQRT( ( ( (R2320 + R2330) * 0.5) / (b*R2230 + a*R2390) ) * ( ( (R2520 + R2530) * 0.5) / (d*R2390 + c*R2600) ) ) )  ;;MISTAKE d was accidently multiplied by 2230 instead of 2390  (CEV 4/12)
     
     
    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return, img

end
