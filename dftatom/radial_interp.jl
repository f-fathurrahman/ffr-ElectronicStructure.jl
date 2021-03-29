# returns the interpolated value of f(x), where f is defined on the
# grid R
# input parameters
# R ... grid
# n ... length of the grid
# f ... function defined on the grid R
# x ... point at which to interpolate
# i ... the integer number of the radial point at which "x" lies
# output parameters
# V ... the interpolated function value
function radial_interp(r, f, x::Float64, i::Int64)
    #IMPLICIT NONE 
    #INTEGER, INTENT(in) :: n
    #REAL(8), INTENT(in) :: f(n), x, R(n)
    #REAL(8), INTENT(out) :: V
    #INTEGER, INTENT(in) :: i
    
    N = size(r,1)
    if N < 4
        error("radial_interp: n >= 4 required")
    end
  
    j1 = i - 1
    j2 = j1 + 1

    n1 = j1 - 1
    n2 = j2 + 1

    if n1 < 1
        n2 = n2 - n1 + 1
        n1 = 1
    end

    if n2 > N
        n1 = n1 - (n2 - N)
        n2 = N
    end
  
    return @views lagrange_interp( r[n1:n2], f[n1:n2], x )
end

