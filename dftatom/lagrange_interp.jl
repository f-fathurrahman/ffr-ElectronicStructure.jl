# Interpolates the x, y points to estimate the value at the point "t".
# val (out) ...... the estimated function value at the point "t".
# t (in).......... the "x" value at which to estimate the value
# X (array, in) .. the x-values of the points
# Y (array, in) .. the y-values of the points
# n (in) ......... the length of the arrays X and Y
# REAL(8), INTENT(in) :: t
# INTEGER, INTENT(in) :: n
# REAL(8), INTENT(in) :: X(n), Y(n)
# REAL(8), INTENT(out) :: val
# REAL(8) :: f, denum
# INTEGER :: i, j
function lagrange_interp(X, Y, xx)
    N = size(X,1)
    retval = 0.0
    for j in 1:N
        f = 1.0
        den = 1.0
        for i in 1:N
            if i == j
                continue
            end
            f = f * ( xx - X[i] )
            den = den * ( X[j] - X[i] )
        end
        retval = retval + Y[j] * f/den
    end
    return retval
end

