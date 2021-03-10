from sympy import *

r = symbols("r")

f = exp(-r/5)
#f = exp(-r*r/5)

s = integrate(f, (r, 0.0, oo))
print("s    = ", s)
print("N(s) = ", N(s))