;------------------------------------------------------------------------
;  Find the 0.64 micron band depth (BD640): (really calculated at 0.648)
;  SLM modified 26 Jun 2007 to move long wavelength shoulder out 
;  of VNIR filter zone boundary, from 680 to 709 nm
;  NOTE: THIS PARAMETER WILL BE DEGRADED BECAUSE THE KEY WAVELENGTH NEAR
;  648 NM IS LOCATION WITHIN THE VNIR FILTER ZONE BOUNDARY
;------------------------------------------------------------------------
function crism_summary_bd640,cube,wvt,hyper=hyper, ignore_val=ignore_val

    return, crism_sumutil_band_depth ( cube, wvt, 600, 648, 709, $ ;
                    hyper=hyper, ignore_val=ignore_val, mid_width=3 )

end
