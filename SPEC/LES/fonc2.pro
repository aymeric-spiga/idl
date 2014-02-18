function fonc2, x, p

y = p(0) * x^(2./3.) * ( 1. - p(1) * x + p(2) )^2
 
y = p(0) * x^(2./3.) * ( p(2) - p(1) * x )^2

y = p(0) * x^(2./3.) * ( 1. - p(1) * x )^2


return, y
end
