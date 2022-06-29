;------------------------------------------------------------------------
;  Find the 2.3 micron Fe-OH band depth (BD2290): (really calculated at 2.293 micron)
;  (with this formulation, this is also a good CO2 ice parameter--2.293
;   micron band depth) (wvs 2290, 2250, 2350)
;------------------------------------------------------------------------
function crism_summary_bd2290,cube,wvt,hyper=hyper,ignore_val=ignore_val

    ; wavelength 2290 is 249 in Zeta2
    ; wavelength 2250 is 255 in Zeta2
    ; wavelength 2350 is 240 in Zeta2

    return, crism_sumutil_band_depth( cube, wvt, 2250, 2290, 2350, $
                    hyper=hyper, ignore_val=ignore_val )     
         
end
