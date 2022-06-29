;11/09/2010 (fps)
;   Initial development - return -1,0,1 for negative, zero, and positive input values
;12/16/2013 (fps)
;   Cleanup for limited release
;   	Don't bother with zero_indx
;   

function sign_ternary, in_thing

;in_thing_size = size(in_thing, /struct)

out_thing = fix(in_thing * 0.0)

negative_indx = where(in_thing LT 0.0, num_negative_indx)
positive_indx = where(in_thing GT 0.0, num_positive_indx)
;zero_indx = where(in_thing EQ 0.0, num_zero_indx)

if (num_negative_indx GT 0) then begin
    out_thing[negative_indx] = -1
endif

if (num_positive_indx GT 0) then begin
    out_thing[positive_indx] = 1
endif

return, out_thing

END
