pro unit_summary_band_names

    ; check to see that object creation fails if no detector value provided
    x = obj_new ( "summary_band_names" )
    if (obj_valid(x)) then begin
        print, "expected object creation to fail without detector argument to constructor"
    endif

    det = 2
    x = obj_new ( "summary_band_names", det )
    if obj_valid (x) then begin

        ; assert type is joint
        if ( x->get_type() ne "joint" ) then begin
            print, "incorrect type value for joint test"
        endif

        ; assert n_elements of names is 46
        names = x->get_names(error=err);
        if ( n_elements(names) ne 46 ) then begin
            print, "incorrect number of names in band name array"
        endif

        ; check that R2700 is at index 42
        y = x->get_index( "R2700", error=err )
        if ( err ne 0 ) then begin
            print, "error returned by get_index()"
        endif else begin
            if ( y ne 42 ) then begin
                print, "band R2700 not at expected index for joint"
            endif
        endelse

        obj_destroy, x
    endif

    ; assert type is long from /ir
    det = 1
    x = obj_new ( "summary_band_names", det )
    if ( x->get_type() ne "long" ) then begin
        print, "incorrect type for long test"
    endif
    if ( x->get_detector() ne det ) then begin
        print, "incorrect detector type for long test"
    endif
    obj_destroy, x

    ; assert type is short from /vnir
    det =  0
    x = obj_new ( "summary_band_names", det )
    if ( x->get_type() ne "short" ) then begin
        print, "incorrect type for short test"
    endif
    if ( x->get_detector() ne det ) then begin
        print, "incorrect detector type for short test"
    endif
    obj_destroy, x



end
