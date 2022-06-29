;
;  Be default, a structure is returned (equivalent to using the /struct switch).
;  If interested in only the resulting cube, use the /data switch
;
;  /extend is used to return not only the summary parameters, but the entire
;  set of bands needed to prepare the browse products.  This has the effect
;  of including an additional set of individual spectral band (or their
;  filtered versions) at the far end of the cube.
;
;10/18/2011 (fps)
;   Modified call to crism_summary_rpeak1() to allow VNIR-only input cubes
;   Papered over RPEAK1-->BDI1000VIS dependency issues - need to revist and veriify for all cases
;12/20/2012 (fps)
;   fixed call to mro_crism_default_wavelengths() to explicitly specify source detector for in_wvt
;01/30/2013 (fps)
;   Commented out VAR parameter call - still needs work...
;05/06/2013 (cev)
;   version=1 is to call the old and new parameters
;   version=0 (default) just calls the new parameters
;   version=3 just calls most 'important' parameters (CEV)
;06/2013 (cev)
;   extend=0 is to not include extended parameters
;   extend=1 (default) calls extended parameters for browse products

function compute_summary_params_cev, in_cube, det, in_wvt, index_vec=index_vec, $
                ignore_val=ignore_val, band_names=band_names, hyper=hyper, $
                extend=extend, data=data, struct=struct, rpeak=rpeak, $
                bdicont=bdicont, version=version

    if ( not (keyword_set(ignore_val))) then begin
        ignore_val = 65535.0
    endif
  
    if n_elements(extend) eq 1 then begin
        extend = extend
    endif else begin
        extend = 1
    endelse
    
    if ( not (keyword_set(version))) then begin
        version = 0
    endif

    t0 = systime(1)

    ; create a fully populated index vector if not provided
    bn = obj_new ( "summary_band_names", det, version, extend=extend)
    band_list = bn->get_names()
    nbands = bn->get_nbands()
    
    if ( not (keyword_set(index_vec))) then begin
        index_vec = replicate(1,nbands)
    endif     
    
    ; return only the band names if so requested
    return_band_names = band_list[where(index_vec eq 1 )]
    if ( keyword_set(band_names)) then begin
        return, return_band_names
    endif
    
    ;check for known bad bands - spectrally subset cube and wavelength vector if necessary
    ;wave_struct = mro_crism_wavelength_defaults(in_wvt, /conservative)
    
    wave_struct = mro_crism_wavelength_defaults(in_wvt, /conservative,ir = (det EQ 1), vnir = (det EQ 0), joined = (det EQ 2))
    wvt = in_wvt[wave_struct.indx]
    cube = in_cube[*,*,wave_struct.indx]
    help, cube
    help, wvt


    ; get set of selected parameter bands from index vector
    selected_indices = where ( index_vec ne 0, nindices )
    if ( nindices gt 0 ) then begin

        ; prepare mask 
        sz = size(cube)
        NB = sz[3]  ; spectral bands
        NX = sz[1]  ; spatial detector
        NY = sz[2]  ; spatial along track

        ; return cube
        pcube = make_array ( nx, ny, nindices, value = ignore_val )

        ; determine list of bands to compute
        selected_bands = band_list[selected_indices]


	;Continuum cals is dependent on BDI1000IR or BDI2000 being selected
	;calc bdicont for input
	if ((max(strpos(selected_bands, 'BDI1000IR')) NE -1) OR (max(strpos(selected_bands, 'BDI2000')) NE -1)) then begin
            bdicont = crism_sumutil_bdicont( cube, wvt, hyper=hyper, ignore_val=ignore_val)
	endif

        ; loop over names of bands selected to be computed
        for iii=0, nindices-1 do begin

            case selected_bands[iii] of

                'R770': begin
                        print, "computing R770"
                        pcube[*,*,iii] = crism_sumutil_single(cube, wvt, 770, hyper=hyper, ignore_val=ignore_val )
                        end
                'RBR': begin
                        print, "computing RBR"
                        pcube[*,*,iii] = crism_summary_rbr(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD530': begin
                        print, "computing BD530"
                        pcube[*,*,iii] = crism_summary_bd530(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD530_2': begin
                        print, "computing BD530_2"
                        pcube[*,*,iii] = crism_summary_bd530_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end        
                'SH600': begin
                        print, "computing SH600"
                        pcube[*,*,iii] = crism_summary_sh600(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'SH600_2': begin
                        print, "computing SH600_2"
                        pcube[*,*,iii] = crism_summary_sh600_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'SH770': begin
                        print, "computing SH770"
                        pcube[*,*,iii] = crism_summary_sh770(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD640': begin
                        print, "computing BD640"
                        pcube[*,*,iii] = crism_summary_bd640(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                
                'BD640_2': begin
                        print, "computing BD640_2"
                        pcube[*,*,iii] = crism_summary_bd640_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD860': begin
                        print, "computing BD860"
                        pcube[*,*,iii] = crism_summary_bd860(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD860_2': begin
                        print, "computing BD860_2"
                        pcube[*,*,iii] = crism_summary_bd860_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'BD920': begin
                        print, "computing BD920"
                        pcube[*,*,iii] = crism_summary_bd920(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD920_2': begin
                        print, "computing BD920_2"
                        pcube[*,*,iii] = crism_summary_bd920_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end			
                'RPEAK1': begin
                        print, "computing RPEAK1"
                        rpeak_struct = crism_summary_rpeak1(cube, wvt, hyper=hyper, ignore_val=ignore_val, joined = (det EQ 2))
                        pcube[*,*,iii] = rpeak_struct.rpeak_wavelength
                        end	
                'BDI1000VIS': begin
                        print, "computing BDI1000VIS"
                        ;pcube[*,*,iii] = crism_summary_bdi1000vis(cube, wvt, hyper=hyper, ignore_val=ignore_val, rpeak=rpeak )
                        pcube[*,*,iii] = crism_summary_bdi1000vis(cube, wvt, hyper=hyper, ignore_val=ignore_val, rpeak=rpeak_struct.rpeak_reflectance )
                        end
                'BDI1000IR': begin
                        print, "computing BDI1000IR"
                        pcube[*,*,iii] = crism_summary_bdi1000ir(cube, wvt, hyper=hyper, ignore_val=ignore_val, bdicont=bdicont )
                        end
                'R1330': begin
                        print, "computing R1330"
                        pcube[*,*,iii] = crism_sumutil_single(cube, wvt, 1330, hyper=hyper, ignore_val=ignore_val, kernel_width=11)
                        end
                'OLINDEX2': begin
                        print, "computing OLINDEX2 "
                        pcube[*,*,iii] = crism_summary_olindex2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'OLINDEX3': begin
                        print, "computing OLINDEX3 "
                        pcube[*,*,iii] = crism_summary_olindex3(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'LCPINDEX': begin
                        print, "computing LCPINDEX"
                        pcube[*,*,iii] = crism_summary_lcpindex(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'LCPINDEX2': begin
                        print, "computing LCPINDEX2"
                        pcube[*,*,iii] = crism_summary_lcpindex2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'HCPINDEX': begin
                        print, "computing HCPINDEX"
                        pcube[*,*,iii] = crism_summary_hcpindex(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'HCPINDEX2': begin
                        print, "computing HCPINDEX2"
                        pcube[*,*,iii] = crism_summary_hcpindex2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end      
                'BD1300': begin
                        print, "computing BD1300"
                        pcube[*,*,iii] = crism_summary_bd1300(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                                             
                'VAR': begin
                        print, "computing VAR"
                        pcube[*,*,iii] = crism_summary_var(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'ISLOPE1': begin
                        print, "computing ISLOPE1"
                        pcube[*,*,iii] = crism_summary_islope1(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD1400': begin
                        print, "computing BD1400"
                        pcube[*,*,iii] = crism_summary_bd1400(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD1435': begin
                        print, "computing BD1435"
                        pcube[*,*,iii] = crism_summary_bd1435(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'BD1500': begin
                        print, "computing BD1500"
                        pcube[*,*,iii] = crism_summary_bd1500(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD1500_2': begin
                        print, "computing BD1500_2"
                        pcube[*,*,iii] = crism_summary_bd1500_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'ICER1': begin
                        print, "computing ICER1" 
                        pcube[*,*,iii] = crism_summary_icer1(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'ICER1_2': begin
                        print, "computing ICER1_2" 
                        pcube[*,*,iii] = crism_summary_icer1_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end        
                'BD1750': begin
                        print, "computing BD1750"
                        pcube[*,*,iii] = crism_summary_bd1750(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD1750_2': begin
                        print, "computing BD1750_2"
                        pcube[*,*,iii] = crism_summary_bd1750_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'BD1900': begin
                        print, "computing BD1900"
                        pcube[*,*,iii] = crism_summary_bd1900(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD1900_2': begin
                        print, "computing BD1900_2"
                        pcube[*,*,iii] = crism_summary_bd1900_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end      
                'BD1900R': begin
                        print, "computing BD1900R"
                        pcube[*,*,iii] = crism_summary_bd1900r(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end    
                'BD1900R2': begin
                        print, "computing BD1900R2"
                        pcube[*,*,iii] = crism_summary_bd1900r2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                                           
                'BDI2000': begin
                        print, "computing BDI2000"
                        pcube[*,*,iii] = crism_summary_bdi2000(cube, wvt, hyper=hyper, ignore_val=ignore_val, bdicont=bdicont )
                        end
                'BD2100': begin
                        print, "computing BD2100"
                        pcube[*,*,iii] = crism_summary_bd2100(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD2100_2': begin
                        print, "computing BD2100_2"
                        pcube[*,*,iii] = crism_summary_bd2100_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end    
                'BD2165': begin
                        print, "computing BD2165"
                        pcube[*,*,iii] = crism_summary_bd2165(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                                             
                'BD2190': begin
                        print, "computing BD2190"
                        pcube[*,*,iii] = crism_summary_bd2190(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end        
                'D2200': begin
                        print, "computing D2200"
                        pcube[*,*,iii] = crism_summary_d2200(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                                                  
                'DOUB2200H': begin
                        ; only compute this parameter if the hyper switch has been set
                        if keyword_set(hyper) then begin
                            print, "computing DOUB2200H"
                            pcube[*,*,iii] = crism_summary_doub2200H(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        endif
                        end      
                'MIN2200': begin
                            print, "computing MIN2200"
                            pcube[*,*,iii] = crism_summary_min2200(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                          
                'BD2210': begin
                        print, "computing BD2210"
                        pcube[*,*,iii] = crism_summary_bd2210(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD2210_2': begin
                        print, "computing BD2210_2"
                        pcube[*,*,iii] = crism_summary_bd2210_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'BD2230': begin
                        print, "computing BD2230"
                        pcube[*,*,iii] = crism_summary_bd2230(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end 
                'MIN2250': begin
                        print, "computing MIN2250"
                        pcube[*,*,iii] = crism_summary_min2250(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end             
                'BD2250': begin
                        print, "computing BD2250"
                        pcube[*,*,iii] = crism_summary_bd2250(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD2265': begin
                            print, "computing BD2265"
                            pcube[*,*,iii] = crism_summary_bd2265(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                       
                'BD2290': begin
                        print, "computing BD2290"
                        pcube[*,*,iii] = crism_summary_bd2290(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'D2300': begin
                        print, "computing D2300"
                        pcube[*,*,iii] = crism_summary_d2300(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD2355': begin
                        print, "computing BD2355"
                        pcube[*,*,iii] = crism_summary_bd2355(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'SINDEX': begin
                        print, "computing SINDEX"
                        pcube[*,*,iii] = crism_summary_sindex(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'SINDEX2': begin
                        print, "computing SINDEX2"
                        pcube[*,*,iii] = crism_summary_sindex2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'ICER2': begin
                        print, "computing ICER2"
                        pcube[*,*,iii] = crism_summary_icer2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'ICER2_2': begin
                        print, "computing ICER2_2"
                        pcube[*,*,iii] = crism_summary_icer2_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'BDCARB': begin
                        print, "computing BDCARB"
                        pcube[*,*,iii] = crism_summary_bdcarb(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'MIN2295_2480': begin
                            print, "computing MIN2295_2480"
                            pcube[*,*,iii] = crism_summary_min2295_2480(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                  
                'MIN2345_2537': begin
                        print, "computing MIN2345_2537"
                        pcube[*,*,iii] = crism_summary_min2345_2537(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                                                
                'BD2500H': begin
                        ; only compute this parameter if the hyper switch has been set
                        if keyword_set(hyper) then begin
                            print, "computing BD2500H"
                            pcube[*,*,iii] = crism_summary_bd2500h(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        endif
                        end
                'BD2500_2': begin
                            print, "computing BD2500_2"
                            pcube[*,*,iii] = crism_summary_bd2500_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'BD3000': begin
                        print, "computing BD3000"
                        pcube[*,*,iii] = crism_summary_bd3000(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD3100': begin
                        print, "computing BD3100"
                        pcube[*,*,iii] = crism_summary_bd3100(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD3200': begin
                        print, "computing BD3200"
                        pcube[*,*,iii] = crism_summary_bd3200(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD3400': begin
                        print, "computing BD3400"
                        pcube[*,*,iii] = crism_summary_bd3400(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'BD3400_2': begin
                        print, "computing BD3400_2"
                        pcube[*,*,iii] = crism_summary_bd3400_2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'CINDEX': begin
                        print, "computing CINDEX"
                        pcube[*,*,iii] = crism_summary_cindex(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'CINDEX2': begin
                        print, "computing CINDEX2"
                        pcube[*,*,iii] = crism_summary_cindex2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end                        
                'BD2600': begin
                        print, "computing BD2600"
                        pcube[*,*,iii] = crism_summary_bd2600(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'R440': begin
                        print, "computing R440"
                        pcube[*,*,iii] = crism_sumutil_single (cube, wvt, 440, hyper=hyper, ignore_val=ignore_val )
                        end
                'IRR1': begin
                        print, "computing IRR1"
                        pcube[*,*,iii] = crism_summary_irr1(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'IRR2': begin
                        print, "computing IRR2"
                        pcube[*,*,iii] = crism_summary_irr2(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                'IRR3': begin
                        print, "computing IRR3"
                        pcube[*,*,iii] = crism_summary_irr3(cube, wvt, hyper=hyper, ignore_val=ignore_val )
                        end
                                
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ;
                ;  All summary bands below this line fall into the category of bands in the extended
                ;  summary parameter cube.  These bands are necessary to populate the browse products,
                ;  but are not official summary parameter bands and only exist in the extended portion
                ;  of the cube
                ;
                'R530': begin
                        print, "computing R530"
                        pcube[*,*,iii] = crism_sumutil_single(cube, wvt, 530, hyper=hyper, ignore_val=ignore_val )
                        end
                'R600': begin
                        print, "computing R600"
                        pcube[*,*,iii] = crism_sumutil_single(cube, wvt, 600, hyper=hyper, ignore_val=ignore_val )
                        end
                'R1080': begin
                        print, "computing R1080"
                        pcube[*,*,iii] = crism_sumutil_single(cube, wvt, 1080, hyper=hyper, ignore_val=ignore_val )
                        end
                'R1506': begin
                        print, "computing R1506"
                        pcube[*,*,iii] = crism_sumutil_single(cube, wvt, 1506, hyper=hyper, ignore_val=ignore_val )
                        end
                'R2529': begin
                        print, "computing R2529"
                        pcube[*,*,iii] = crism_sumutil_single(cube, wvt, 2529, hyper=hyper, ignore_val=ignore_val )
                        end
                'R3920': begin
                        print, "computing R3920"
                        pcube[*,*,iii] = crism_sumutil_single(cube, wvt, 3920, hyper=hyper, ignore_val=ignore_val )
                        end
                else: stop
            endcase
        endfor
    endif

    print, "Time Spent computing summary parameters = ", systime(1)-t0
    help, pcube

    if ( keyword_set(data) ) then begin
        return, pcube
    endif else begin
        return, { sumparams:pcube, band_names:return_band_names }
    endelse
    
end
