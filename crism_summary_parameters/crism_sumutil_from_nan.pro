;
;  This is a simple function that will replace any IEEE NAN values
;  (ie:  !values.f_nan) with the value provided as an input variable
;  (ie:  CRISM_NAN or 65535) 
;
function crism_sumutil_from_nan, input, newval

    output = input
    w = where (finite( input, /nan), count)
print, "CRISM_SUMUTIL_FROM_NAN:  #replaced = "+strtrim(count,2)
    if (count gt 0)  then begin
        output[w] = newval
    endif
    return, output

end
