function lanczos(pw, Potentials, k)
#
# Usage: [T,V,f] = lanczos(H,v0,k)
#
# Purpose:
#    Perform k-step Lanczos iterations
#
# Input:
#    H  (Ham object)       Hamiltonian
#    v0 (Wavefun object)   Wavefunction
#    k  (integer)          The number of Lanczos iterations
#
# Output:
#  T  (k by k matrix)    The projection of H in 
#                        the Krylov subspace spanned by columns of V;
#  V  (Wavefun object)   Stores the Lanczos vectors;
#  f  (Wavefun object)   Stores the residual vector;
#
#  This version of Lanzos performs full (re)orthogonalization
#  instead of using a 3-term recurrence. It doesn't check
#  for early convergence (i.e., an invariant subspace is found
#  before the kth step of Lanczos is completed.)
#
    Ngwx = pw.gvec.Ngwx
    v = rand(Ngwx,1)
    #
    T = zeros(k,k)
    #
    beta = norm(v)
    v = v ./ beta
    #
    V = zeros(Ngwx,k)
    V[:,1] = v
    Hv = op_H( pw, Potentials, v )
    h = V[:,1]' * Hv
    #
    T[1,1] = h
    f = Hv - V[:,1]*h
    # One-step of reorthogonalization
    s = V[:,1]' * f
    h = h + s
    f = f - V[:,1]*s
    # MAIN LOOP
    for j = 2:k
        beta = norm(f)
        T(j,j-1) = beta
        v = f/beta
        V = [V v]
        Hv = op_H(pw, Potentials, v)
        h = V' * Hv
        f = Hv - V*h
        # One-step of reorthogonalization
        s = V' * f
        h = h + s
        f = f - V*s
        #
        T(1:j,j) = h
    end

end