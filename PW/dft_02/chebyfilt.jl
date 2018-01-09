function  chebyfilt(pw, Potentials, X, degree, lb, ub)
    Ngwx    = size(X)[1]
    Nstates = size(X)[2]
    #
    e = (ub - lb)/2
    c = (ub + lb)/2
    sigma = e/(lb-ub)
    sigma1 = sigma
    #
    Y = zeros(Complex128,Ngwx,Nstates)
    Y1 = zeros(Complex128,Ngwx,Nstates)
    #
    Y = op_H(pw, Potentials, X) - X*c
    Y = Y*sigma1/e
    #
    for i = 2:degree
        sigma2 = 1/(2/sigma1 - sigma)
        Y1 = ( op_H(pw,Potentials,Y) - Y*c)*2 * sigma2/e - X*(sigma*sigma2)
        X = Y
        Y = Y1
        sigma = sigma2
    end
    return Y
end