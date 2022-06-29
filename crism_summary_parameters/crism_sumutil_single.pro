;
;  This routine extracts a single band from a cube.  If the hyper switch
;  is not provided, it simply extracts the requested band.  If the hyper switch
;  is set, if computes an image that is median filtered in only the spectral
;  direction using the neighboring "kernel width" bands.
;
;  To ensure that the median and smooth functions work properly, any crism_nan
;  values are first replaced with the IDL IEEE NAN values.
;
;  If the /norestore switch is provided, IEEE NAN values are left in the return
;  image.  By default, IEEE NAN values are replaced with the provided
;  ignore_val value (typically CRISM_NAN = 65535).
;
;  
;  Howard Taylor
;  JHU/APL
;  Aug 8, 2011
;
function crism_sumutil_single,cube,wvt,band_wavelength, hyper=hyper, $
        ignore_val=ignore_val, kernel_width=kernel_width, avg=avg,  $
        norestore=norestore

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; default to a kernel of length 5 in the spectral direction
    if (n_elements(kernel_width) eq 0) then begin
        kernel_width= 5
    endif

    ; identify index into cube for specific wavelength
    R_indx=mro_crism_lookupwv(band_wavelength,wvt)

    ; Filter the data if hyper switch set, else choose individual band
    if keyword_set(hyper) then begin

        ; identify adjacent indices matching +/- half a kernel width from the base wavelength
        half_width = fix(kernel_width / 2)
        minidx = ( R_indx - half_width ) > 0  ; clamp to 0 if R-half_width is less than 0
        maxidx = ( R_indx + half_width ) < (n_elements(wvt)-1)  ; clamp to max index if R+half_width is past the end
        

        ; check the wavelength extremes to ensure there are no gaps in wavelength within the kernel
        wavelength_coverage = abs( wvt[minidx] - wvt[maxidx])
        mean_wavelength_per_channel = 6.55 ; nm/channel
        max_allowable_missing_channels = 2.0 
        wavelength_gap_tolerance = mean_wavelength_per_channel * $
                        ( kernel_width + max_allowable_missing_channels )
        if ( wavelength_coverage gt wavelength_gap_tolerance ) then begin
            print, "kernel wavelength coverage exceeds max wavelength gap tolerance. ( " + $
                   strtrim(wavelength_coverage,2)+" nm < " + $
                   strtrim(wavelength_gap_tolerance,2)+" nm )"
            print, "    kernel width = "+strtrim(kernel_width,2)
            print, "    center wavelength = "+strtrim(wvt[R_indx],2)
            print, "    wavelength table index = "+strtrim(R_indx,2)
            fmtstr="(A,"+strtrim((kernel_width-1),2)+"(F8.2,','),F8.2,A)"
            print, "    kernel wavelengths: [", wvt[minidx:maxidx], "]", format=fmtstr
        endif

        ; create a subsetted wavelength table and find the corresponding index in the subsetted table
        wvt_n=wvt[minidx:maxidx]  
            
        R_indx_n=(mro_crism_lookupwv(band_wavelength,wvt_n))[0]

        ; subset the "kernel_width" spectral pixels centered on the base wavelength and
        ; replace CRISM_NAN with IEEE floating point NAN
        subcube = crism_sumutil_to_nan( cube[*,*,minidx:maxidx], ignore_val)
        
        ; filter the subsetted cube using either boxcar or median filtering
        if ( keyword_set(avg)) then begin
            smthcube = smooth( subcube, [1,1,kernel_width], /edge_truncate, /nan)
            img = smthcube[*,*,R_indx_n]
        endif else begin
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; CEV edits
            if half_width eq 0 then begin
            img=subcube
            endif else begin
            img = median(subcube, dimension=3, /even)
            endelse
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        endelse

    endif else begin

        ; extract a single band from the input cube and replace CRISM_NAN with IEEE floating point NAN
        img = reform(crism_sumutil_to_nan( cube[*,*,R_indx], ignore_val))

    endelse

    ; By default, replace any lingering IEEE NaN values with CRISM_NAN.  If
    ;  norestore is set, leave the IEEE NaN values in the return img
    if not(keyword_set(norestore)) then begin
        img = crism_sumutil_from_nan(img, ignore_val )
    endif

    return, img
end
