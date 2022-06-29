;------------------------------------------------------------------------
;  Find the 0.86 micron band depth (BD860): ('hematite band')
;  SLM modified 26 Jun 2007 to move long wavelength shoulder 
;  from 920 to 984 nm to sample more spectral curvature
;------------------------------------------------------------------------
function crism_summary_bd860,cube,wvt,hyper=hyper,ignore_val=ignore_val

    return, crism_sumutil_band_depth( cube, wvt, 800, 860, 984, $  ;
                    hyper=hyper, ignore_val=ignore_val, mid_width=11 )

end
