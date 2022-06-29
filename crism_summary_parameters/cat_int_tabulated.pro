function cat_int_tabulated, x, f

; IDL int_tabulated with filtering of non-unique x values. 
; Used by summary parameter integrations to guard against the 
; case where bad bands might cause parameter wavelengths to 
; resolve to the same band.

s = sort(x)
xs = x[s]
u = s[uniq(xs)]

; Inputs with no repeated X values:
xx = x[u]
ff = f[u]

return, int_tabulated(xx,ff)
end
