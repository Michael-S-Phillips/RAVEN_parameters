;------------------------------------------------------------------------
;Find spectral variance parameter (VAR):
; SLM modified 26 Jun 2007 to remove 1630 and 1660 nm which may have 
; significant artifacts due to their location at the IR zone 1 zone 2 
; filter boundary
; CEV cleaned up to put hyperspectral formulation back in and
; make sure CRISM Nan is maintained (5/13/2013)
;------------------------------------------------------------------------

function crism_summary_var,cube,wvt,hyper=hyper,ignore_val=ignore_val

    if (not(keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif

    sz = size(cube)
    NX = sz[1]  ; spatial detector ( 640 detector columns unless binned )
    NY = sz[2]  ; spatial # of frames along track
    NB = sz[3]  ; spectral bands ( # wavelengths )

    ;From Omega:
    ;Fit a line from 1-2.3 microns and find variance of observed values from
    ; fit values by summing in quadrature over 26 intervening wavelengths.
    ; (note: ol & px will have high variance b/c of bands--other stuff should
    ; have low variance.)
    ;Find reflectances that we want to fit a line to:
    ;For OMEGA: VARSLOPELAMS=['0.983450','2.286900']
    ; IN CRISM, the shortest IR Row is 445 with a wavelength of around 1001.6 nm
    ; in Zeta2 wavelength table, so we use that
    ; Closest to 2.286 microns  is row 250, 2285.95 nm
    ; 
    ; Get intervening wavelengths and observed reflectances at those values:
    ; 26 intervening wavelengths in Omega are
    ; ;LAMS=['1.026400','1.055000','1.083700','1.126800','1.213000','1.328000',$
    ;        '1.328000','1.371100','1.399900','1.428600','1.471700','1.471700',$
    ;        '1.514800','1.700700','1.814300','1.814300','1.856800','1.927200',$
    ;        '1.983400','2.011300','2.067000','2.122500','2.150100','2.163900',$
    ;        '2.205100','2.232400','2.246100']
    ; In the CRISM wavetables:
    
    if keyword_set(hyper) then begin
        R1021_indx = mro_crism_lookupwv ( 1021, wvt )
        R2253_indx = mro_crism_lookupwv ( 2253, wvt )
        indices = [ R1021_indx, R2253_indx]
        wavelength_indx = indgen (max(indices) - min(indices) + 1) + min(indices)
        wvs=wvt[wavelength_indx]

    endif else begin
        wvs_in = [ 1021, 1048, 1080, 1153, 1212, 1258, 1278, 1330, 1370, 1396, 1429, $
                1469, 1508, 1561, 1693, 1752, 1812, 1878, 1930, 1983, 2009, 2068, $
                2121, 2141, 2167, 2206, 2233, 2253]

        wavelength_indx = make_array( n_elements(wvs_in), value = ignore_val)
    
        for i = 0, n_elements(wvs_in)-1 do begin
          wavelength_indx[i]=mro_crism_lookupwv(wvs_in[i],wvt)
         
        endfor

    endelse
    
    wvs=wvt[wavelength_indx]
     ; extract the related wavelengths from the cube, replacing CRISM NaN with IEEE NaN
    var_cube = crism_sumutil_to_nan( cube, ignore_val)
    
    varslopelams = [mro_crism_lookupwv(1014,wvt,/w), $
                    mro_crism_lookupwv(2287,wvt,/w) ]

    varslopewl   = [mro_crism_lookupwv(1014,wvt), $
                    mro_crism_lookupwv(2287,wvt)]
                    
    varslopers = fltarr ( nx, ny, 2 ) ; Not sure why there is a 2 here...(CEV)

    for k = 0 , n_elements(varslopelams)-1 do begin
        varslopers[*,*,k] = float( var_cube [ *, *, varslopewl[k]] )
    endfor    
 
    ; Find the actual reflectances
    obsrs = make_array( nx, ny, n_elements(wvs), value = ignore_val)
    for  k = 0, n_elements(wvs)-1 do begin
        indx = mro_crism_lookupwv( wvs[k], wvt)
        obsrs[ *, *, k] = float( var_cube[*,*,indx])
    endfor

    ;  Fit a line and find the variance:
    var = make_array( nx, ny, value=ignore_val)
    for j= 0, ny-1 do begin
        for i=0, nx-1 do begin
            fit = linfit( float(varslopelams), varslopers[i,j,*] )
            predrs = fit[0] + fit[1]*wvs
            var[i,j] = total( (predrs - obsrs[i,j,*])^2 )
        endfor
    endfor
    
        ; replace the IEEE NAN values with CRISM_NAN
     var = crism_sumutil_from_nan( var, ignore_val)

    return, var

end

