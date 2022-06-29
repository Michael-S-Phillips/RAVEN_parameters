;  Find the 2.35 micron CO band depth (BD2350): (wvs 2320, 2330, 2350, 2290, 2430)
function crism_summary_bd2350,cube,wvt,hyper=hyper,ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; get single bands from the cube
    R2290 = crism_sumutil_single(cube, wvt, 2290, hyper=hyper, /norestore ); 228 in Zeta2
    R2320 = crism_sumutil_single(cube, wvt, 2320, hyper=hyper, /norestore )
    R2330 = crism_sumutil_single(cube, wvt, 2330, hyper=hyper, /norestore )
    R2350 = crism_sumutil_single(cube, wvt, 2350, hyper=hyper, /norestore ); 249 in Zeta2
    R2430 = crism_sumutil_single(cube, wvt, 2430, hyper=hyper, /norestore )

    ; get wavelengths of those closest to those requested
    W2290 = (mro_crism_lookupwv( 2290, wvt ))[0]
    W2320 = (mro_crism_lookupwv( 2320, wvt ))[0]
    W2330 = (mro_crism_lookupwv( 2330, wvt ))[0]
    W2350 = (mro_crism_lookupwv( 2350, wvt ))[0]
    W2430 = (mro_crism_lookupwv( 2430, wvt ))[0]

    WL = W2290
    WC = ( W2290 + W2330 + W2350 ) / 3.0
    WH = W2430
    a = 0.33333
    b = 0.33333
    c = 0.33333
    d = (WC-WL)/(WH-WL)     ; d gets multipled by the longer (higher wvln)  band
    ee = 1.0-d                ; ee gets multiplied by the shorter (lower wvln) band

    ; compute bd2350 using IEEE NaNs
    img = 1.0 - ((a*R2320 + b*R2330 + c*R2350)/ (ee*R2290 + d*R2430))

    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img
end
