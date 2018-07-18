pairs(n::Int64) = ((i, j) for i = 1:n for j = 1:i)
rpairs(n::Int64) = ((i,j) for i in 1:n for j in 1:n) # rectangular option to old pairs
spairs(n::Int64) = ((i, j) for i = 1:n for j = 1:(i-1)) # subdiagonal option to old pairs
 
triangle(i::Int64) = div(i*(i+1),2)
triangle(i::Int64,j::Int64) = i<j ? triangle(j-1)+i : triangle(i-1)+j
                        
iiterator(n::Int64) = ((i,j,k,l) for (i,j) in pairs(n) for (k,l) in pairs(n) if triangle(i,j) <= triangle(k,l))

iindex(i::Int64,j::Int64,k::Int64,l::Int64) = triangle(triangle(i,j),triangle(k,l))
trace2(A,B) = sum(A.*B)

#=
function all_1e_ints( bfs::BasisSet, atoms::Atoms )
    n = length(bfs)
    S = Array{Float64}(undef,n,n)
    T = Array{Float64}(undef,n,n)
    V = Array{Float64}(undef,n,n)
    for j = 1:n
        for i = j:n
            S[i,j] = overlap(bfs[i], bfs[j])
            S[j,i] = S[i,j]
            T[i,j] = kinetic(bfs[i], bfs[j])
            T[j,i] = T[i,j]
            V[i,j] = nuclear_attraction(bfs[i], bfs[j], atoms)
            V[j,i] = V[i,j]
        end
    end
    return S,T,V
end
=#

function all_1e_ints(bfs::BasisSet, atoms::Atoms)
    n = length(bfs)
    S = Array{Float64}(undef,n,n)
    T = Array{Float64}(undef,n,n)
    V = Array{Float64}(undef,n,n)
    for (i,j) in pairs(n)
        a,b = bfs[i], bfs[j]
        S[i,j] = S[j,i] = overlap(a,b)
        T[i,j] = T[j,i] = kinetic(a,b)
        V[i,j] = V[j,i] = nuclear_attraction(a,b,atoms)
    end
    return S,T,V
end

#=
function all_twoe_ints(bfs,ERI=coulomb)
    n = length(bfs)
    totlen = div(n*(n+1)*(n*n+n+2),8)
    ints2e = zeros(totlen)  # using zeros instead of using Array constructor
    for i = 1:n
    for j = 1:n
        ij = triangle(i-1,j)
        for k = 1:n
        for l = 1:n
            kl = triangle(k-1,l)
            if ij <= kl
                ints2e[iindex(i,j,k,l)] = ERI(bfs[i], bfs[j], bfs[k], bfs[l])
            end
        end
        end
    end
    end
    return ints2e
end
=#

function all_twoe_ints(bfs,ERI=coulomb)
    n = length(bfs)
    totlen = div(n*(n+1)*(n*n+n+2),8)
    ints2e = Array{Float64}(undef,totlen)
    for (i,j,k,l) in iiterator(n)
        ints2e[iindex(i,j,k,l)] = ERI(bfs[i], bfs[j], bfs[k], bfs[l])
    end
    return ints2e
end

#=
function make2JmK(D::Array{Float64,2},Ints::Array{Float64,1})
    n = size(D,1)
    G = Array{Float64}(undef,n,n)
    D1 = reshape(D,n*n)
    temp = Array{Float64}(undef,n*n)
    for j = 1:n
    for i = j:n
        kl = 1
        for k = 1:n
            start = 1
            for l = start:n
                temp[kl] = 2*Ints[iindex(i,j,k,l)]-Ints[iindex(i,k,j,l)]
                kl += 1
            end
        end
        G[i,j] = dot(D1,temp)
        G[j,i] = G[i,j]
    end
    end
    return G
end
=#

function make2JmK(D::Array{Float64,2},Ints::Array{Float64,1})
    n = size(D,1)
    G = Array{Float64}(undef,n,n)
    D1 = reshape(D,n*n)
    temp = Array{Float64}(undef,n*n)
    for (i,j) in pairs(n)
        kl = 1
        for (k,l) in rpairs(n)
            temp[kl] = 2*Ints[iindex(i,j,k,l)]-Ints[iindex(i,k,j,l)]
            kl += 1
        end
        G[i,j] = G[j,i] = dot(D1,temp)
    end
    return G
end

dmat(U::Array{Float64,2},nocc::Int64) = U[:,1:nocc]*U[:,1:nocc]'

function RHF( atoms::Atoms, MaxIter::Int64=20, verbose::Bool=true)
    @printf("Starting RHF calculation:")

    @printf("Building basis ...")
    bfs = build_basis(atoms)
    @printf("... done\n\n")

    @printf("1e integrals ...")
    S,T,V = all_1e_ints(bfs,atoms)
    @printf("... done\n\n")

    @printf("2e integrals ...")
    Ints = all_twoe_ints(bfs)
    @printf("... done\n\n")

    h = T+V
    E,U = eigen(h,S)
    println("E = ", E)
    Enuke = nuclear_repulsion(atoms)
    Nelectrons = get_Nelectrons(atoms)
    nclosed,nopen = divrem( Int64(Nelectrons), 2 )
    Eold = 0.0
    Energy = 0.0

    @printf("Nel = %d, Nclosed = %d, nopen = %d\n", Nelectrons, nclosed, nopen)
  
    if verbose
        println("S=\n$S")
        println("h=\n$h")
        println("T=\n$T")
        println("V=\n$V")
        println("E: $E")
        println("U: $U")
        println("2e ints:\n$Ints")
    end

    IS_CONVERGED = false
    dEtot = 1.0
    for iter in 1:MaxIter
        D = dmat(U,nclosed)
        #if verbose
        #    println("D=\n$D")
        #end
        G = make2JmK(D,Ints)
        H = h + G
        E,U = eigen(H,S)
        Eone = trace2(D,h)
        Etwo = trace2(D,H)
        Energy = Enuke + Eone + Etwo
        dEtot = abs(Energy - Eold)
        #println("HF: $iter  $Energy : $Enuke  $Eone  $Etwo $dEtot")
        @printf("HF: %4d %18.10f %18.10f %18.10f: %18.10f %18.10e\n",
                iter, Enuke, Eone, Etwo, Energy, dEtot)
        if isapprox(Energy,Eold)
            IS_CONVERGED = true
            break
        end
        Eold  = Energy
    end

    if !IS_CONVERGED
        println("WARNING: RHF is not converged")
    end

    return Energy,E,U
end
