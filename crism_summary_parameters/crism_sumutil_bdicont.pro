;------------------------------------------------------------------------
;
;  This function models the continuum from 1um to 2um and is used for both
;  the bdi1000ir and bdi2000 calculations.
;
;  Old comments from mro_crism_summary_params.pro:
;
;  03/19/2008    (fps) 
;    Total rewrite of the BDI1000IR calculation
;      Revised continuum calculation - used for both BDI1000IR and BDI2000
;      Find median value at multispectral wavelengths over the interval ~[1330,1810]
;      Find median value at multispectral wavelengths over the interval ~[2430,2600]
;      Linear fit between the median (wavelength, value) points is the continuum
;   Revised ~1000 nm integrated band depth calculation
;       Integration of normalized reflectance values over the inteveral ~[1045,1255]
;  04/02/2008    (fps)
;    Refined continuum calculation used in BDI1000IR and BDI2000 calculations
;      Find upper quartile (75th percentile) value at multispectral wavelengths over the interval ~[1330,1810]
;      Linear fit between the ~2500 nm median (wavelength, value) and ~1500 nm upper quartile (wavelength, value) points
;  05/09/2013 (cev)
;    Changed wavelength fit for ~1500 nm to a median (rather than upper quartile).   
;------------------------------------------------------------------------

function crism_sumutil_bdicont, cube, wvt, hyper=hyper, ignore_val=ignore_val


    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; determine cube size for later use
    sz = size(cube)
    nx = sz[1]  ; x = spatial (640 detector columns, unless binned)
    ny = sz[2]  ; y = spatial (# of frames along track)
    nb = sz[3]  ; z = spectral # wavelengths 

    ; Identify candidate wavelengths to model spectrum near 1 um for continuum calculation
    if keyword_set(hyper) then begin
        R1330_indx = mro_crism_lookupwv(1330,wvt)
        R1875_indx = mro_crism_lookupwv(1875,wvt)
        indices = [R1330_indx, R1875_indx]
        wavelength_indx = indgen (max(indices) - min(indices) + 1) + min(indices)

    endif else begin
        ; pick off wavelengths from multispectral sampling
        R1330_indx = mro_crism_lookupwv(1330,wvt)
        R1370_indx = mro_crism_lookupwv(1370,wvt)
        R1395_indx = mro_crism_lookupwv(1395,wvt)
        R1425_indx = mro_crism_lookupwv(1425,wvt)
        R1465_indx = mro_crism_lookupwv(1465,wvt)
        R1500_indx = mro_crism_lookupwv(1500,wvt)
        R1505_indx = mro_crism_lookupwv(1505,wvt)
        R1560_indx = mro_crism_lookupwv(1560,wvt)
        R1625_indx = mro_crism_lookupwv(1625,wvt)
        R1660_indx = mro_crism_lookupwv(1660,wvt)
        R1690_indx = mro_crism_lookupwv(1690,wvt)
        R1750_indx = mro_crism_lookupwv(1750,wvt)
        R1810_indx = mro_crism_lookupwv(1810,wvt)
        R1875_indx = mro_crism_lookupwv(1875,wvt)
        wavelength_indx = [R1330_indx,R1370_indx,R1395_indx,R1425_indx,$
                           R1465_indx,R1500_indx,R1505_indx,R1560_indx,$
                           R1625_indx,R1660_indx,R1690_indx,R1750_indx,$
                           R1810_indx]


    endelse

    wavelength_vec = wvt[wavelength_indx]
    ;print, '  Wavelength vector:'
    ;print, wavelength_vec
    ;print, '  Index vector:'
    ;print, wavelength_indx

    ; subset the cube, replacing all CRISM NaN values with IEEE Nan values
    cont1_cube = crism_sumutil_to_nan( cube[*,*,wavelength_indx], ignore_val)

    ; compute the percentile statistics from the subsetted cube
    percentile1_frame = make_array(nx, ny, value = ignore_val)         
    percentile1_wavelength = make_array(nx,ny, value = ignore_val)
    ;help, percentile1_frame            
    ;help, percentile1_wavelength

    ; compute upper quartile frame and corresponding wavelength frame
    for i = 0, (nx - 1) do begin
        for j = 0, (ny - 1) do begin
        
            spec_vec = cont1_cube[i,j,*]

            ; IEEE NaNs may be in the data here.  No longer checking for CRISM NaNs
            sort_indx = sort(spec_vec)
            select_indx = sort_indx[n_elements(wavelength_indx) * 0.75] ;~upper quartile indx
            
            percentile1_frame[i,j] = spec_vec[select_indx]
            ;percentile1_wavelength[i,j] = wavelength_vec[select_indx]  ;~upper quartile wvt
            percentile1_wavelength[i,j]=median(wavelength_vec, /even)  ;median wvt -- changed by CEV (05/09/2013)
            
        endfor
    endfor

    ;median1_frame = median(cont1_cube, dimension = 3)
    ;help, median1_frame
    ;median1_wavelength = median(wavelength_vec, /even)
    ;help, median1_wavelength

    ; compute spectral response near 2.4 um for long response of continuum
    if keyword_set (hyper) then begin
        R2430_indx = mro_crism_lookupwv(2430,wvt)
        R2600_indx = mro_crism_lookupwv(2600,wvt)
        indices = [R2430_indx, R2600_indx]
        wavelength_indx = indgen (max(indices) - min(indices) + 1) + min(indices)
    endif else begin
        R2430_indx = mro_crism_lookupwv(2430,wvt)
        R2460_indx = mro_crism_lookupwv(2460,wvt)
        R2530_indx = mro_crism_lookupwv(2530,wvt)
        R2600_indx = mro_crism_lookupwv(2600,wvt)
        wavelength_indx = [R2430_indx,R2460_indx,R2530_indx,R2600_indx]
    endelse
    wavelength_vec = wvt[wavelength_indx]
                        
    ; subset the cube for the end of the continuum near 2.4 um
    ;  replacing CRISM_NaNs with IEEE NaNs
    cont2_cube = crism_sumutil_to_nan( cube[*,*,wavelength_indx], ignore_val)
    median2_frame = median(cont2_cube, dimension = 3)
    median2_wavelength = median(wavelength_vec, /even)
    cont_slope_frame = (median2_frame - percentile1_frame) / $
                            (median2_wavelength - percentile1_wavelength)

    ; return a structure with enough details to represent the continuum
    return, { cont_slope_frame:cont_slope_frame,        $
              median2_frame:median2_frame,              $
              median2_wavelength:median2_wavelength,    $
              percentile1_frame:percentile1_frame,      $
              percentile1_wavelength:percentile1_wavelength }

end
