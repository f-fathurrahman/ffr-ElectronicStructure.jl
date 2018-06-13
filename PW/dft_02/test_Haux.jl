include("../common/ortho_gram_schmidt.jl")
include("printMatrix.jl")

function test_main()
    srand(2222)
    Ngwx = 1000
    Nstates = 4
    psi = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)
    psi_ortho1 = ortho_gram_schmidt(psi)

    U = inv(sqrtm(psi'*psi))
    printMatrix(U)

    U2 = (psi*U)'*(psi*U)
    printMatrix(U2)

    U1 = psi_ortho1'*psi_ortho1
    printMatrix(U1)

    Haux = zeros(ComplexF64,Nstates,Nstates)

end

test_main()