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
;7/1/12 (cev)
;Utilize most available multispectral wavelengths to 'integrate' under the fit
;Fit shifted slightly to remove contributions from HCP to the OLINDEX
;Because this is the only parameter that extrapolates the estimated 'continuum' (linear fit) past
;the to tie points for the fit, it has the possibility of creating negative 'Rc' values.  Therefore,
;divide by the absolute value of Rc, rather than just Rc, in order to preserve the negative. (Values below zero
;hold no meaning, since the index will approach 0 and the Rc value becomes more negative).
;
;3/27/13 (cev)
;Scaled weighting factors on calculation so that they add to 1.
;
;11/17/13 (cevb)
;One last major overhaul before delivery.  Changed so that the anchors are at shorter wavelengths (1860 and 1750 nm).
;Makes OLINDEX3 distinctly *less sensitive to hydrated minerals that tend to have a drop off past 2.3 microns.
;False positives include:  CO2 ice and sometimes gypsum and jarosite.  OLINDEX3 best detects forsterite (Mg-rich olivine)
;while BD1300 best detects fayalite (Fe-rich olivine).

function crism_summary_olindex3, cube, wvt, hyper=hyper, ignore_val=ignore_val, classic = classic

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif 


    ; extract channels from cube replacing CRISM_NAN with IEEE_NaN

    R1210 = crism_sumutil_single(cube, wvt, 1210, hyper=hyper, /norestore, kernel_width = 7) 
    R1250 = crism_sumutil_single(cube, wvt, 1250, hyper=hyper, /norestore, kernel_width = 7)
    R1263 = crism_sumutil_single(cube, wvt, 1263, hyper=hyper, /norestore, kernel_width = 7)
    R1276 = crism_sumutil_single(cube, wvt, 1276, hyper=hyper, /norestore, kernel_width = 7)
    R1330 = crism_sumutil_single(cube, wvt, 1330, hyper=hyper, /norestore, kernel_width = 7) 
    
    R1750 = crism_sumutil_single(cube, wvt, 1750, hyper=hyper, /norestore, kernel_width = 7) 
    R1862 = crism_sumutil_single(cube, wvt, 1862, hyper=hyper, /norestore, kernel_width = 7) 
    
    ; identify nearest CRISM wavelength

    W1210 = (mro_crism_lookupwv(1210,wvt,/w))[0]
    W1250 = (mro_crism_lookupwv(1250,wvt,/w))[0]
    W1263 = (mro_crism_lookupwv(1263,wvt,/w))[0]
    W1276 = (mro_crism_lookupwv(1276,wvt,/w))[0]
    W1330 = (mro_crism_lookupwv(1330,wvt,/w))[0]

    W1750 = (mro_crism_lookupwv(1750,wvt,/w))[0]
    W1862 = (mro_crism_lookupwv(1862,wvt,/w))[0] 
    
    ; compute the corrected reflectance interpolating 
    slope = ( R1862 - R1750 ) / ( W1862 - W1750 )   ;slope = ( R2120 - R1690 ) / ( W2120 - W1690 )
    intercept = R1862 - slope * W1862               ;intercept = R2120 - slope * W2120

    ; weighted sum of relative differences

    Rc1210 = slope * W1210 + intercept
    Rc1250 = slope * W1250 + intercept
    Rc1263 = slope * W1263 + intercept
    Rc1276 = slope * W1276 + intercept
    Rc1330 = slope * W1330 + intercept

     img = (((Rc1210-R1210)/(abs(Rc1210))) * 0.1) + (((Rc1250-R1250)/(abs(Rc1250))) * 0.1) + $
      (((Rc1263-R1263)/(abs(Rc1263))) * 0.2) + (((Rc1276-R1276)/(abs(Rc1276))) * 0.2) + (((Rc1330-R1330)/(abs(Rc1330))) * 0.4)   
               
               
    ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan ( img, ignore_val)

    return,img

end
