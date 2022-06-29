;------------------------------------------------------------------------
; VER3 Addition: Carbonate index (CINDEX):
; (used as a replacement for BD3800 because CRISM won't go out to 4050nm
; CRISM ROW 1 (3925 nm) is very noisy in groundcal. Subsitute row 2
; (3919.38 nm) for better result for now. (wvs 3750, 3630, 3950)
;------------------------------------------------------------------------
function crism_summary_cindex,cube,wvt, hyper=hyper, ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract channels from cube replacing CRISM_NAN with IEEE_NaN
    R3630 = crism_sumutil_single(cube, wvt, 3630, $     ; add width keyword
                    hyper=hyper, /norestore)  ; 46 in Zeta2
    R3750 = crism_sumutil_single(cube, wvt, 3750,  $    ; add width keyword
                    hyper=hyper, /norestore)  ; 28 in Zeta2
    R3920 = crism_sumutil_single(cube, wvt, 3920,  $    ; add width keyword
                    hyper=hyper, /norestore)  ;  2 in Zeta2

    ; identify nearest CRISM wavelength
    W3630 = (mro_crism_lookupwv(3630,wvt,/w))[0]
    W3750 = (mro_crism_lookupwv(3750,wvt,/w))[0]
    W3920 = (mro_crism_lookupwv(3920,wvt,/w))[0]

    ;  compute the one-sided band depth calculation for the radiance value at
    ;  3920 by estimating the continuum using the slope given by the change in
    ;  radiance between R3630 and R3750.  Use this change in radiance, along with
    ;  the change in wavelength to estimate a slope for the continuum.  Given the
    ;  radiance at R3750, extrapolate the continuum value at R3920 given the
    ;  change in wavelength from 3750 to 3920 and the slope computed previously to
    ;  estimate the continuum.  Then measure the one-sided band depth as the ratio
    ;  at a wavelength of 3920 compared to the continuum estimate there.
    ;  
    img = (R3750 + ((R3750 - R3630)/(W3750 - W3630)) * (W3920 - W3750)) / R3920 - 1
    
    ;;checked above and calculated correctly

    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img

end
