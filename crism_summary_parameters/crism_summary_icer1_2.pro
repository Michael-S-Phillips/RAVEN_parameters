function crism_summary_icer1_2,cube,wvt,hyper=hyper,ignore_val=ignore_val
                
           img1= crism_summary_bd1435(cube, wvt, hyper=hyper, ignore_val=ignore_val )
           img2= crism_summary_bd1500_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )

    img=(img2 ne 65535.)*(1.-((1.-img1)/(1.-img2)))+(img2 eq 65535.)*65535.
    
  ; replace the IEEE NAN values with CRISM_NAN
    img = crism_sumutil_from_nan(img, ignore_val)
    
        return, img

end
