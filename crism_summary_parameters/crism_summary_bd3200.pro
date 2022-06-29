;------------------------------------------------------------------------
;  Find the 3.2 micron H2O ice band depth (BD3200): (really calculaed at 3.32 micron)
;   (wvs 3320, 3250, 3390)
;------------------------------------------------------------------------
function crism_summary_bd3200,cube,wvt,hyper=hyper,ignore_val=ignore_val

    ; wavelength 3250 is 104  in Zeta2
    ; wavelength 3320 is  93 in Zeta2
    ; wavelength 3390 is  82 in Zeta2

    return, crism_sumutil_band_depth( cube, wvt, 3250, 3320, 3390, $
                    hyper=hyper, ignore_val=ignore_val )

end
