;
;  This is a simple function that will replace any occurrences of
;  the specified value (ie CRISM_NAN or 65535 ) replace it with 
;  the IEEE NaN value (ie:  !values.f_nan) to simplify calculations
;  to prevent the need for exhaustive scans for CRISM_NAN values
;
function crism_sumutil_to_nan, input, crism_nan

    output = input
    w = where ( input eq crism_nan, count )
    print, "CRISM_SUMUTIL_TO_NAN:  #replaced = "+strtrim(count,2)
    if (count gt 0)  then begin
        ; even if input is a double array, this will be promoted from float
        output[w] = !values.f_nan
    endif
    return, output

end
