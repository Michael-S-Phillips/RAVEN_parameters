;
;  summary parameter cube band names object
;
;  case-insensitive search for band_names
;
;  Howard Taylor
;  JHU/APL
;  2011 May 10
;
;
;05/07/2013 (cev & fps)  
;   added a new method get_version() for the parameter list version
;   added a new tag in the band names object structure for the version
;06/05/2013 (cev)
;   Changed 'extend' call to match with compute_summary_params_cev
;   so that extend=1 adds extra parameters for browse products
;   Changed ordering of parameters by increasing wvl of feature
;12/05/2013 (fps)
;   Added detector discrimination to /extend band name set
;


;======================================================
;  return list of names as an array
;======================================================
function summary_band_names::get_names, error=err

    names='unknown'
    ;
    ;  Check to make sure that names is available
    ;
    if ( ptr_valid(self.names) eq 0 ) then begin
        err=1
        goto, get_names_skip
    endif else err=0

    names = make_array( self.nbands, /string, value='')
    for iii=0,(self.nbands-1) do begin
        names[iii] = (*self.names)[iii]
    endfor

get_names_skip:
    return, names

end


;======================================================
;  return name of a specific index
;======================================================
function summary_band_names::get_name, indx, error=err

    name='unknown'
    ;
    ;  Check to make sure that names is available
    ;
    if ( ptr_valid(self.names) eq 0 ) then begin
        err=1
        goto, get_name_skip
    endif else err=0

    ; check that indx is in bounds
    if (( indx lt 0 ) or (indx gt (self.nbands-1))) then begin
        goto, get_name_skip
    endif

    ; pick off the name at the specified index
    name = (*self.names)[indx]

get_name_skip:
    return, name

end

;======================================================
;  return index of a specific band name
;======================================================
function summary_band_names::get_index, search_name, error=err

    indx=-1
    ;
    ;  Check to make sure that names is available
    ;
    if ( ptr_valid(self.names) eq 0 ) then begin
        err=1
        goto, get_index_skip
    endif else err=0

    ; find the search_name in the list of names
    indx = where (strlowcase(*self.names) eq strlowcase(search_name), count )
    if ( count eq  0) then err = 1  else err = 0

get_index_skip:
    return, indx

end

;======================================================
;  return type for the cube (short, long, joint)
;======================================================
function summary_band_names::get_type
    return, self.type
end

;======================================================
;  return nbands for the cube
;======================================================
function summary_band_names::get_nbands
    return, self.nbands
end

function summary_band_names::get_detector
    return, self.detector
end

function summary_band_names::get_version 
    return, self.version
end
   
;======================================================
;  Natural ordering of short summary parameter band names
;======================================================
function summary_band_names::short_summary_band_names
    ;names = [   'R770', 'RBR', 'BD530', 'SH600', 'SH600_2', 'SH770', $
    ;            'BD640', 'BD860', 'BD860_2', 'BD920', 'BD920_2', $
    ;            'RPEAK1', 'BDI1000VIS', 'R440', 'IRR1' ] 
    
    v=self->get_version()
    case v of
    1: begin
    names = [   'R770', 'RBR', 'BD530', 'BD530_2', 'SH600', 'SH600_2', $
                'SH770', 'BD640', 'BD640_2', 'BD860', 'BD860_2', 'BD920', $
                'BD920_2', 'RPEAK1', 'BDI1000VIS', 'R440', 'IRR1' ] 
                ; CEV May/12, added SH600_2, BD920_2, BD860_2
                ; CEV July/12, added SH770
                ; CEV Feb/13 added BD530_2 - now for real!
                ; CEV March/13 added BD640_2
    end
    0: begin
        names = [   'R770', 'RBR', 'BD530_2', 'SH600_2', $
                'SH770', 'BD640_2', 'BD860_2', $
                'BD920_2', 'RPEAK1', 'BDI1000VIS', 'R440', 'IRR1' ] 
    end
    
    3: begin
        names = [   'BD530_2',  $
                'BD640_2', 'BD860_2', $
                'BD920_2' ] 
    end
    endcase     
    return, names
end

;======================================================
;  Natural ordering of long summary parameter band names
;======================================================
function summary_band_names::long_summary_band_names
 ;   names = [   'BDI1000IR', 'IRA', 'OLINDEX2', 'OLINDEX3', 'LCPINDEX', 'LCPINDEX2', 'HCPINDEX', 'HCPINDEX2', 'VAR',  $
 ;               'ISLOPE1', 'BD1435', 'BD1500', 'BD1500_2', 'ICER1', 'BD1750', 'BD1750_2', 'BD1900', 'BD1900_2',    $
 ;               'BDI2000', 'BD2100', 'BD2210', 'BD2100_2', 'BD2190', 'OPALINDEX', 'DOUB2200H', 'BD2265H', 'BD2290', 'D2300', 'SINDEX', 'SINDEX2',     $
 ;               'ICER2', 'BDCARB', 'BDCARB2', 'BDCARBMGH', 'BD3000', 'BD3100', 'BD3200', 'BD3400', 'BD3400_2',     $
 ;               'CINDEX', 'CINDEX2', 'BD1270O2', 'BD1400H2O', 'BD2000CO2', 'BD2350',       $
 ;               'BD2600', 'IRR2', 'R2700', 'BD2700', 'IRR3', 'BD2500H', 'BD2500H2'  ]  
    
    v=self->get_version()
    case v of
    1: begin
    names = [   'BDI1000IR', 'R1330', 'OLINDEX2', 'OLINDEX3', 'LCPINDEX', 'LCPINDEX2', $
                'HCPINDEX', 'HCPINDEX2', 'BD1300', 'VAR', 'ISLOPE1', 'BD1400', 'BD1435', 'BD1500', $
                'BD1500_2', 'ICER1', 'ICER1_2', 'BD1750', 'BD1750_2', 'BD1900', 'BD1900_2', 'BD1900R', 'BD1900R2',  $
                'BDI2000', 'BD2100', 'BD2100_2', 'BD2165', 'BD2190', 'DOUB2200H', 'MIN2200', $
                'BD2210', 'BD2210_2', 'D2200', 'BD2230','BD2250', 'MIN2250', 'BD2265', 'BD2290', 'D2300', 'BD2355', 'SINDEX', $
                'SINDEX2', 'ICER2', 'ICER2_2', 'BDCARB', 'MIN2295_2480', 'MIN2345_2537', 'BD2500H', $
                'BD2500_2', 'BD3000', 'BD3100', 'BD3200', 'BD3400', 'BD3400_2', $
                'CINDEX', 'CINDEX2', 'BD2600', 'IRR2', 'IRR3' ]
                ; CEV May/12, added LCPINDEX2, HCPINDEX2, BD1500_2, BD1750_2, BD1900_2, BD2100_2, BD2190, OPALINDEX, DOUB2200H, BD2250, BD2265H, 
                ; SINDEX2, BDCARB2, BDCARBMGH, BD3400_2, CINDEX2, and BD2500H2
                ; CEV Jan/13, added BD1300, BD1400, BD2355, changed OPALINDEX to DOUB2250
                ; CEV 4/1/13 added BD1900R, BD1900R2, BD2165, DOUB2200H2
                ; CEV 4/10/13 DOUB2200H2->MIN2200H, DOUB2250->MIN2250, BDCARB2->MIN2342and2537, BDCARBMGH->MIN2302and2505
                ; CEV 5/8/13 added BD2230
                ; CEV 5/15/13 renamed MIN2200H->MIN2200, MIN2347and2537->MIN2345_2537, MIN2302and2505->MIN2295_2480, BD2500H2->BD2500_2, BD2265H->BD2265
                ; CEV 5/21/13 added D2200
       end
    0: begin
        names = [   'BDI1000IR', 'R1330', 'OLINDEX3', 'LCPINDEX2', $
                'HCPINDEX2', 'BD1300', 'VAR', 'ISLOPE1', 'BD1400', 'BD1435', $
                'BD1500_2', 'ICER1_2', 'BD1750_2', 'BD1900_2', 'BD1900R2',  $
                'BDI2000', 'BD2100_2', 'BD2165', 'BD2190', 'MIN2200', $
                'BD2210_2', 'D2200', 'BD2230', 'BD2250', 'MIN2250', 'BD2265', 'BD2290', 'D2300', 'BD2355', $
                'SINDEX2', 'ICER2_2', 'MIN2295_2480', 'MIN2345_2537',  $
                'BD2500_2', 'BD3000', 'BD3100', 'BD3200', 'BD3400_2', $
                'CINDEX2', 'BD2600', 'IRR2', 'IRR3' ]
       end
       
     3: begin
        names = [  'OLINDEX3', 'LCPINDEX2', $
                'HCPINDEX2', 'BD1300', 'BD1400', 'BD1435', $
                'BD1500_2','BD1750_2', 'BD1900R2',  $
                'BD2100_2', 'BD2165', 'BD2190', 'MIN2200', $
                'BD2210_2', 'D2200', 'BD2230', 'BD2250', 'MIN2250', 'BD2265', 'BD2290', 'D2300', 'BD2355', $
                'SINDEX2', 'MIN2295_2480', 'MIN2345_2537', $
                'BD2500_2' ]
       end
       
    endcase
    return, names
end

;======================================================
;  Natural ordering of joint summary parameter band names 
;======================================================
function summary_band_names::joint_summary_band_names
        names = [ self->short_summary_band_names(), self->long_summary_band_names() ]
    return, names
end

function summary_band_names::joint_extend, det
    ;print, det
    case (det) of
	(0):begin   ;vnir
    	    names = [ 'R530', 'R600' ]
	end

	(1):begin   ;ir
       	    names = [ 'R1080', 'R1506', 'R2529', 'R3920' ]
	end

	else:begin  ;joined or indeterminate
            names = [ 'R530', 'R600', 'R1080', 'R1506', 'R2529', 'R3920' ]
	endelse
    endcase

    return, names
end

;======================================================
;  object constructor
;======================================================
function summary_band_names::init, det, version, extend=extend

    if (n_elements(det) eq 0 ) then goto, init_errors
    if (n_elements(version) eq 0 ) then goto, init_errors
    
    self.version = version
    self.detector = det
      
    ;  Create a mapping from index to band name
    case det of
        0: begin
            self.type = "short" ; vnir
            list = self->short_summary_band_names()
            end
        1: begin
            self.type = "long" ; ir
            list = self->long_summary_band_names()
            end
        2: begin
            self.type = "joint"
            list = self->joint_summary_band_names()
            end
        else: goto, init_errors
    endcase

    ;
    ;  If the /extend switch is specified, tack on other bands required for
    ;  browse products
    ;
    if extend eq 1 then begin
        list = [ list, self->joint_extend(det)]
    endif
   
    self.names = ptr_new(list)
    self.nbands = n_elements(list)

    return, 1   ; success

init_errors:
    return, 0   ; failure

end

;======================================================
;  object deconstructor
;======================================================
pro summary_band_names::cleanup
    if ( ptr_valid (self.names )) then ptr_free, self.names
end

;======================================================
;  band names object structure
;======================================================
pro summary_band_names__define

    tmp = { summary_band_names, $
            names: ptr_new(),   $
            nbands: 0,          $
            detector: 0,        $
            version: 0,         $
            type: ''            $
    }
end
