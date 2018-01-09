function  Y = chebyfilt(H, X, degree, lb, ub)
%
% Y = chebyfilt(H, X, degree, lb, ub);
%
% Purpose:
%   Apply p(H) to X to obtain Y = p(H)*X, where p(lambda) is a shifted
%   and scaled Chebyshev polynomial of degree 'degree'.  The polynomial
%   is bounded between -1 and 1 within [lb,ub].
%
% Input:
%   H (Ham object)      Hamiltonian 
%   X (Wavefun object)  Wavefunction
%   degree (integer)    The degree of the Chebyshev polynomial
%   lb (float)          The lower bound of the interval in which the 
%                       polynomial is bounded between -1 and 1.
%   ub (float)          The upper bound of the interval in which the 
%                       polynomial is bounded between -1 and 1.
%
% Output:  
%
%   Y  (Wavefun object)  Wavefunction Y = p(H)*X;
%
e = (ub-lb)/2;
c = (ub+lb)/2;
sigma = e/(lb-ub);
sigma1 = sigma;
Y = H*X - X*c;
Y = Y*sigma1/e;
for i = 2:degree
   sigma2 = 1/(2/sigma1 - sigma);
   Y1 = (H*Y-Y*c)*2*sigma2/e-X*(sigma*sigma2);
   X = Y;
   Y = Y1;
   sigma = sigma2;
end;
