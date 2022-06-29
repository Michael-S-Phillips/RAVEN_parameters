;------------------------------------------------------------------------
;
;  This function is used to compute generic band depth
;
;  A continuum is determined from the DN value of the cube at the low
;  and high wavelengths.  The parameter is then computed by examining
;  the ratio between the value at the center wavelength compared to the
;  interpolated value of the continuum at that center wavelength.
;
;  If the /hyper switch is specified, each of these DN values is then
;  computed by taking other neighboring spectral values into consideration.
;  This sampling of the center wavelength is either done by a simple average
;  of its neighbors or by using a median of a neighboring set.
;
;  The size of the population of the kennel width can be provided, or
;  a default value is used otherwise.  ( see crism_sumutil_single()).
;
;  For each of the low, middle or high wavelengths, specificy these:
;
;       low_width:  width of the kernel (spectrally) for average or median
;       low_avg:    use boxcar averaging rather than default of median
;
;  Howard Taylor
;  JHU/APL
;  06/03/2011
;
;------------------------------------------------------------------------
function crism_sumutil_band_depth,cube,wvt,low,mid,hi,hyper=hyper, $
            low_width=low_width, low_avg=low_avg, $
            mid_width=mid_width, mid_avg=mid_avg, $
            hi_width=hi_width,   hi_avg=hi_avg,   $
            ignore_val=ignore_val
               
    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract single bands from the cube, replacing CRISM_NAN with IEEE NAN
    Rlow = crism_sumutil_single(cube, wvt, low, hyper=hyper, /norestore, $
        ignore_val=ignore_val, avg=low_avg, kernel_width=low_width )
    Rmid = crism_sumutil_single(cube, wvt, mid, hyper=hyper, /norestore, $
        ignore_val=ignore_val, avg=mid_avg, kernel_width=mid_width )
    Rhi = crism_sumutil_single(cube, wvt, hi, hyper=hyper, /norestore, $
        ignore_val=ignore_val, avg=hi_avg, kernel_width=hi_width )

    ; determine wavelength values for closest crism channels
    WL = (mro_crism_lookupwv(low,wvt,/w))[0]
    WC = (mro_crism_lookupwv(mid,wvt,/w))[0]
    WH = (mro_crism_lookupwv(hi, wvt,/w))[0]
    a = (WC-WL)/(WH-WL)     ; a gets multipled by the longer band
    b = 1.0-a               ; b gets multiplied by the shorter band

    ; compute the band depth using precomputed a and b
    img = 1.0 - ( Rmid / ( b * Rlow + a * Rhi))

    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan( img, ignore_val)

    return,img

end
