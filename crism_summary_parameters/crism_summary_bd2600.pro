;------------------------------------------------------------------------
;  Find the 2.6 micron H2O band depth (BD2600): VER3 Formulation (wvs 2600, 2530, 2630)
;------------------------------------------------------------------------
function crism_summary_bd2600,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 2530, 2600, 2630, $
                    hyper=hyper, ignore_val=ignore_val )

end
