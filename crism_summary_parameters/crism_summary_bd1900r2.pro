;------------------------------------------------------------------------
; Find the 1.90 micron H2O band depth (BD1900R): (wvs 1908-1941 (in), 1862-1875, 2112-2126)
; (CEV) Continuum removed ratios so it is not longer sensitive to slope effects
;------------------------------------------------------------------------

function crism_summary_bd1900r2,cube,wvt,hyper=hyper,ignore_val=ignore_val

if keyword_set(hyper) then begin

      ; extract individual channels, replacing CRISM_NANs with IEEE_NaNs
    R1908 = crism_sumutil_single(cube, wvt, 1908, hyper=hyper, /norestore, kernel_width = 1) 
    R1914 = crism_sumutil_single(cube, wvt, 1914, hyper=hyper, /norestore, kernel_width = 1) 
    R1921 = crism_sumutil_single(cube, wvt, 1921, hyper=hyper, /norestore, kernel_width = 1) 
    R1928 = crism_sumutil_single(cube, wvt, 1928, hyper=hyper, /norestore, kernel_width = 1) 
    R1934 = crism_sumutil_single(cube, wvt, 1934, hyper=hyper, /norestore, kernel_width = 1) 
    R1941 = crism_sumutil_single(cube, wvt, 1941, hyper=hyper, /norestore, kernel_width = 1) 
    
    R1862 = crism_sumutil_single(cube, wvt, 1862, hyper=hyper, /norestore, kernel_width = 1) 
    R1869 = crism_sumutil_single(cube, wvt, 1869, hyper=hyper, /norestore, kernel_width = 1) 
    R1875 = crism_sumutil_single(cube, wvt, 1875, hyper=hyper, /norestore, kernel_width = 1) 
    R2112 = crism_sumutil_single(cube, wvt, 2112, hyper=hyper, /norestore, kernel_width = 1) 
    R2120 = crism_sumutil_single(cube, wvt, 2120, hyper=hyper, /norestore, kernel_width = 1) 
    R2126 = crism_sumutil_single(cube, wvt, 2126, hyper=hyper, /norestore, kernel_width = 1) 
    
    R1815 = crism_sumutil_single(cube, wvt, 1815, hyper=hyper, /norestore);
    R2132 = crism_sumutil_single(cube, wvt, 2132, hyper=hyper, /norestore); 
    
        ; retrieve the CRISM wavelengths nearest the requested values
    W1908 = (mro_crism_lookupwv ( 1908, wvt, /w ))[0] 
    W1914 = (mro_crism_lookupwv ( 1914, wvt, /w ))[0]
    W1921 = (mro_crism_lookupwv ( 1921, wvt, /w ))[0] 
    W1928 = (mro_crism_lookupwv ( 1928, wvt, /w ))[0] 
    W1934 = (mro_crism_lookupwv ( 1934, wvt, /w ))[0] 
    W1941 = (mro_crism_lookupwv ( 1941, wvt, /w ))[0]
    
    W1862 = (mro_crism_lookupwv ( 1862, wvt, /w ))[0]
    W1869 = (mro_crism_lookupwv ( 1869, wvt, /w ))[0]
    W1875 = (mro_crism_lookupwv ( 1875, wvt, /w ))[0]
    W2112 = (mro_crism_lookupwv ( 2112, wvt, /w ))[0]
    W2120 = (mro_crism_lookupwv ( 2120, wvt, /w ))[0]
    W2126 = (mro_crism_lookupwv ( 2126, wvt, /w ))[0] 
    
    W1815 = (mro_crism_lookupwv ( 1815, wvt, /w ))[0]    
    W2132 = (mro_crism_lookupwv ( 2132, wvt, /w ))[0]   
    
        ; compute the interpolated continuum values at selected wavelengths between 1815 and 2530
    slope = ( R2132 - R1815 ) / ( W2132 - W1815 )
    CR1908 = R1815 + slope * ( W1908 - W1815 )
    CR1914 = R1815 + slope * ( W1914 - W1815 )
    CR1921 = R1815 + slope * ( W1921 - W1815 ) 
    CR1928 = R1815 + slope * ( W1928 - W1815 )
    CR1934 = R1815 + slope * ( W1934 - W1815 )
    CR1941 = R1815 + slope * ( W1941 - W1815 ) 
       
    CR1862 = R1815 + slope * ( W1862 - W1815 )
    CR1869 = R1815 + slope * ( W1869 - W1815 )
    CR1875 = R1815 + slope * ( W1875 - W1815 )    
    CR2112 = R1815 + slope * ( W2112 - W1815 )
    CR2120 = R1815 + slope * ( W2120 - W1815 )
    CR2126 = R1815 + slope * ( W2126 - W1815 )
    
    img= 1.0-( (R1908/CR1908+R1914/CR1914+R1921/CR1921+R1928/CR1928+R1934/CR1934+R1941/CR1941) / $
            (R1862/CR1862+R1869/CR1869+R1875/CR1875+R2112/CR2112+R2120/CR2120+R2126/CR2126) )
endif else begin
    R1928 = crism_sumutil_single(cube, wvt, 1928, hyper=hyper, /norestore, kernel_width = 1)
    R1875 = crism_sumutil_single(cube, wvt, 1875, hyper=hyper, /norestore, kernel_width = 1)
    R2112 = crism_sumutil_single(cube, wvt, 2112, hyper=hyper, /norestore, kernel_width = 1) 
    
    R1815 = crism_sumutil_single(cube, wvt, 1815, hyper=hyper, /norestore);
    R2132 = crism_sumutil_single(cube, wvt, 2132, hyper=hyper, /norestore);    
    
    W1928 = (mro_crism_lookupwv ( 1928, wvt, /w ))[0]
    W1875 = (mro_crism_lookupwv ( 1875, wvt, /w ))[0]
    W2112 = (mro_crism_lookupwv ( 2112, wvt, /w ))[0]  
      
    W1815 = (mro_crism_lookupwv ( 1815, wvt, /w ))[0]    
    W2132 = (mro_crism_lookupwv ( 2132, wvt, /w ))[0]    
     
    slope = ( R2132 - R1815 ) / ( W2132 - W1815 )    
    
    CR1928 = R1815 + slope * ( W1928 - W1815 )    
    CR1875 = R1815 + slope * ( W1875 - W1815 )    
    CR2112 = R1815 + slope * ( W2112 - W1815 )    
    
    img= 1.0-( 2.*(R1928/CR1928) / ((R1875/CR1875)+(R2112/CR2112)) )    
    
endelse 
  
  ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan(img, ignore_val)

return, img

end
