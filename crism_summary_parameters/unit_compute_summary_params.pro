function unit_compute_summary_params

    ir_image = "TESTING/FRT0000955B_07_RA164L_TRR3.IMG"
    ir_label = "TESTING/FRT0000955B_07_RA164L_TRR3.LBL"
    vnir_image = "TESTING/FRT0000955B_07_RA164S_TRR3.IMG"
    vnir_label = "TESTING/FRT0000955B_07_RA164S_TRR3.LBL"
    joint_image = "TESTING/FRT0000ABB5_07_IF164J_TRR3.IMG"      ; MSL site with mafics, primary and altered mineralogy
    joint_label = "TESTING/FRT0000ABB5_07_IF164J_TRR3.LBS"

    joint_image = "TESTING/FRT000081B1_07_IF164J_TRR3.IMG"      ; GLO example
    joint_label = "TESTING/FRT000081B1_07_IF164J_TRR3.LBS"

    joint_image = "TESTING/FRT000052BC_07_IF163J_TRR3.IMG"      ; ICE2 example
    joint_label = "TESTING/FRT000052BC_07_IF163J_TRR3.LBS"

    joint_image = "TESTING/FRT000064D9_07_IF166J_TRR3.IMG"      ; TRU, VNA, FEM, FM2, IRA, FAL example  ( good for D2300 too)
    joint_label = "TESTING/FRT000064D9_07_IF166J_TRR3.LBS"

    joint_image = "TESTING/FRT000089F7_07_IF166J_TRR3.IMG"      ; bad band at the top?
    joint_label = "TESTING/FRT000089F7_07_IF166J_TRR3.LBS"

    joint_image = "TESTING/FRT000064D9_07_IF166J_TRR3.IMG"      ; good for OLINDEX2
    joint_label = "TESTING/FRT000064D9_07_IF166J_TRR3.LBS"

    joint_image = "TESTING/FRT000084FA_07_IF163J_TRR3.IMG"      ; good for HCP
    joint_label = "TESTING/FRT000084FA_07_IF163J_TRR3.LBS"
    
    joint_image = "TESTING/FRT00008F86_07_IF163J_TRR3.IMG"      ; good for HCP
    joint_label = "TESTING/FRT00008F86_07_IF163J_TRR3.LBS"

    ;cube = mro_crism_quick_read(vnir_image, /load_struct)
    ;detector = 0

    ;cube = mro_crism_quick_read(ir_image, /load_struct)
    ;detector = 1

    cube = mro_crism_quick_read(joint_image, /load_struct)
    detector = 2

    band_names = compute_summary_params ( cube.image, detector, cube.wavelength_vector, /band_names, /extend )
    ;print, band_names

    ; create an index vector with length matching band_names
    index_vec=  replicate(0, n_elements(band_names))

    index_vec[where ( band_names eq "R770" )] = 1
    index_vec[where ( band_names eq "RBR" )] = 1
    index_vec[where ( band_names eq "R440" )] = 1
    index_vec[where ( band_names eq "BD530" )] = 1
    index_vec[where ( band_names eq "BD920" )] = 1
    index_vec[where ( band_names eq "BD1435" )] = 1     ; significant striping
    index_vec[where ( band_names eq "BD1500" )] = 1     ; checked
    index_vec[where ( band_names eq "BD1900" )] = 1
    index_vec[where ( band_names eq "BD2100" )] = 1     ; some striping
    index_vec[where ( band_names eq "BD2210" )] = 1    ; significant striping
    index_vec[where ( band_names eq "SINDEX" )] = 1    ; faint striping
    index_vec[where ( band_names eq "HCPINDEX" )] = 1
    index_vec[where ( band_names eq "LCPINDEX" )] = 1
    index_vec[where ( band_names eq "OLINDEX2" )] = 1
    index_vec[ where ( band_names eq "D2300" ) ] = 1   ; faint striping
    index_vec[ where ( band_names eq "IRA" ) ] = 1
    index_vec[ where ( band_names eq "SH600" ) ] = 1
    index_vec[ where ( band_names eq "IRR1" ) ] = 1     ; verify with Frank on boxcar params for hyper (max=3e-5)
    index_vec[ where ( band_names eq "IRR2" ) ] = 1     ; verify with Frank on boxcar params for hyper
    index_vec[ where ( band_names eq "IRR3" ) ] = 1     ; verify with Frank on boxcar params for hyper ; stripy (1.1-1.4)
    index_vec[ where ( band_names eq "ICER1" ) ] = 1     ; verify with Frank on boxcar params for hyper ; stripy (1.1-1.4)
    index_vec[ where ( band_names eq "ICER2" ) ] = 1     ; verify with Frank on boxcar params for hyper ; stripy (1.1-1.4)
    index_vec[ where ( band_names eq "R530" ) ] = 1
    index_vec[ where ( band_names eq "R600" ) ] = 1
    index_vec[ where ( band_names eq "R1080" ) ] = 1
    index_vec[ where ( band_names eq "R1506" ) ] = 1
    index_vec[ where ( band_names eq "R2529" ) ] = 1
    index_vec[ where ( band_names eq "R2700" ) ] = 1
    index_vec[ where ( band_names eq "R3920" ) ] = 1
    index_vec[ where ( band_names eq "ISLOPE1" ) ] = 1     ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD2350" ) ] = 1      ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD3000" ) ] = 1      ; checked
    index_vec[ where ( band_names eq "RPEAK1" ) ] = 1      ; verify implementation with Frank
    index_vec[ where ( band_names eq "BDI1000VIS" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD640" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD860" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD920" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD1750" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD2290" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD3100" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD3200" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD3400" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD2600" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD2700" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BDCARB" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BDI1000IR" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BDI2000" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "CINDEX" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "VAR" ) ] = 1  ; verify implementation with Frank
    index_vec[ where ( band_names eq "BD2500H" ) ] = 1  ; verify implementation with Frank

    hyper=1
    extend=1
    x = compute_summary_params ( cube.image, detector, cube.wavelength_vector, $
                    index_vec=index_vec, hyper=hyper, extend=extend, /struct )
    print, "Band Names for returned Summary Parameter Cube:"
    print, x.band_names
    look, x.sumparams

    return, x
end
