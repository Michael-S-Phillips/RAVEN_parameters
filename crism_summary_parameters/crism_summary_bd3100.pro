;------------------------------------------------------------------------
;  Find the 3.1 micron H2O ice band depth (BD3100): (wvs 3120, 3000, 3250)
;------------------------------------------------------------------------
function crism_summary_bd3100,cube,wvt,hyper=hyper,ignore_val=ignore_val

    ; wavelength 3000 is ??? in Zeta2
    ; wavelength 3120 is 123 in Zeta2
    ; wavelength 3250 is 104 in Zeta2

    return, crism_sumutil_band_depth( cube, wvt, 3000, 3120, 3250, $
                    hyper=hyper, ignore_val=ignore_val )

end
