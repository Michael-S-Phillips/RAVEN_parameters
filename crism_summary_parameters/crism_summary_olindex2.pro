;09/30/2011 (fps)
;   Working to resolve OLINDEX2 parameter distribution strangeness...
;   Reference for OLINDEX2 is Salvatore et. al (2010) JGR E07005. It appears that:
;   	1) The description (section 2 - Methods) of the OLINDEX2 equation (Table 2) in the paper is not completely accurate
;   	2) The original MTRDR pipeline implementation was accurate w/r/t the description, but not the equation.
;   The core issue is whether the *difference* or *relative difference* between the expected and observed reflectance values at the 
;   	OLINDEX reference wavelengths are multiplied by the canonical coefficients and the result totaled. The description implies a 
;   	weighted sum of differences, while the equation implies a weighted sum of relative differences.
;   Note: There is an extra open parenthesis in the OLINDEX2 equation in the paper
;   Note: The OLINDEX equation in the paper is incomplete - missing distribution zero-recentering (-1) term
;
;   Also need to resolve 1080 nm vs. 1054 nm disconnect between OLINDEX and OLINDEX2
;
;   Original MTRDR implementation is accessible with keyword /classic
;

function crism_summary_olindex2, cube, wvt, hyper=hyper, ignore_val=ignore_val, classic = classic

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    ; extract channels from cube replacing CRISM_NAN with IEEE_NaN
    R1080 = crism_sumutil_single(cube, wvt, 1080, hyper=hyper, /norestore, kernel_width = 7) 
    R1210 = crism_sumutil_single(cube, wvt, 1210, hyper=hyper, /norestore, kernel_width = 7) 
    R1330 = crism_sumutil_single(cube, wvt, 1330, hyper=hyper, /norestore, kernel_width = 7) 
    R1470 = crism_sumutil_single(cube, wvt, 1470, hyper=hyper, /norestore, kernel_width = 7) 
    
    R1750 = crism_sumutil_single(cube, wvt, 1750, hyper=hyper, /norestore, kernel_width = 7) 
    R2400 = crism_sumutil_single(cube, wvt, 2400, hyper=hyper, /norestore, kernel_width = 7) 

    ; identify nearest CRISM wavelength
    W1080 = (mro_crism_lookupwv(1080,wvt,/w))[0]
    W1210 = (mro_crism_lookupwv(1210,wvt,/w))[0]
    W1330 = (mro_crism_lookupwv(1330,wvt,/w))[0]
    W1470 = (mro_crism_lookupwv(1470,wvt,/w))[0]
    
    W1750 = (mro_crism_lookupwv(1750,wvt,/w))[0]  ;2120? 1750?
    W2400 = (mro_crism_lookupwv(2400,wvt,/w))[0]

    ; compute the corrected reflectance interpolating 
    slope = ( R2400 - R1750 ) / ( W2400 - W1750 )
    intercept = R2400 - slope * W2400

if (keyword_set(classic)) then begin
; ***** Original MTRDR implementation - weighted sum of differences
    bbd1080 = slope * W1080 + intercept - R1080
    bbd1210 = slope * W1210 + intercept - R1210
    bbd1330 = slope * W1330 + intercept - R1330
    bbd1470 = slope * W1470 + intercept - R1470

    img = bbd1080*0.1 + bbd1210*0.1 + bbd1330*0.4 + bbd1470*0.4

endif else begin
; ***** Reference equation implementation - weighted sum of relative differences
    Rc1080 = slope * W1080 + intercept
    Rc1210 = slope * W1210 + intercept
    Rc1330 = slope * W1330 + intercept
    Rc1470 = slope * W1470 + intercept

    ;img = (1-(R1080/Rc1080) * 0.10) + (1-(R1210/Rc1210) * 0.10) + (1-(R1330/Rc1330) * 0.40) + (1-(R1470/Rc1470) * 0.40)
    ; order of operations requires extra parentheses (CEV 3/27/13):
    
    img = ((1-(R1080/Rc1080)) * 0.10) + ((1-(R1210/Rc1210)) * 0.10) + ((1-(R1330/Rc1330)) * 0.40) + ((1-(R1470/Rc1470)) * 0.40)  
        
endelse

    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img

end
