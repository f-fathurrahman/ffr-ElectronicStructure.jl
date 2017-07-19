import PyPlot
const plt = PyPlot

function eval_LF1d_p( α::Int64, L::Float64, grid_x::Array{Float64}, x)
    Nbasis = size(grid_x)[1]
    N = (Nbasis-1)/2
    f = 0.0
    for l = 1:Nbasis
        k = -N + (l-1)
        f = f + sqrt(1/L/Nbasis)*cos(2*pi*k*(x-grid_x[α])/L)
    end
    return f
end

function gen_T_matrix( L::Float64, grid_x::Array{Float64} )
    Nbasis = size(grid_x)[1]
    N = (Nbasis-1)/2
    # Transformation matrix
    T = Array{Complex128}(Nbasis,Nbasis)
    for α = 1:Nbasis
        for l = 1:Nbasis
            k = -N + (l-1)
            x = grid_x[α]
            T[α,l] = sqrt(1/Nbasis)*exp(2*pi*im*k*x/L)
        end
    end
    return T
end

# Example periodic function
function myfunc(L::Float64, x::Float64)
    ω = 2*pi/L
    f = cos(ω*x)*sin(2*ω*x)
    return f
end

function do_plot_myfunc( NptsPlot::Int64, L::Float64, lim )
    NptsPlot = 200
    x = Array{Float64}(linspace(lim[1], lim[2],NptsPlot));
    y = Array{Float64}(NptsPlot)
    for i = 1:NptsPlot
        y[i] = myfunc(L,x[i])
    end
    plt.clf()
    plt.grid()
    plt.plot(x,y)
    plt.savefig("myfunc.png", dpi=300)
end

function eval_from_ex_coefs( ex_coef::Array{Float64},
        L::Float64, grid_x::Array{Float64}, x::Float64 )
    Nbasis = size(ex_coef)[1]
    f = 0.0
    for i = 1:Nbasis
        f = f + ex_coef[i]*eval_LF1d_p( i, L, grid_x, x )
    end
    return f
end

function eval_FBR( l::Int64, L::Float64, N::Int64, x )
    k = -N + (l-1)
    f = 1/sqrt(L) * exp(im*2*pi*k*x/L)
    return f
end

function FBR_eval_from_ex_coefs( ex_coef, L, x )
    Nbasis = size(ex_coef)[1]
    N = Int64( (Nbasis-1)/2 )
    f = 0.0 + im*0.0
    for i = 1:Nbasis
        f = f + ex_coef[i]*eval_FBR( i, L, N, x )
    end
    return f
end

#function eval_DVR_from_FBR( α, T, x )
#    Nbasis = size(T)[1]
#    f = 0.0 + im*0.0
#    for l = 1:Nbasis
#        f = f + T'[l,α]*eval_FBR()
#    end
#end

function main(L::Float64, N::Int64)
    @printf("L = %f\n", L)
    @printf("N = %d\n", N)
    Nbasis = 2*N + 1
    @printf("Nbasis = %d\n", Nbasis)
    # Create grid
    grid_x = Array{Float64}(Nbasis)
    for α = 1:Nbasis
        grid_x[α] = L/Nbasis*(α-N-1)
    end
    Δ = L/Nbasis

    f1 = eval_LF1d_p( 1, L, grid_x, grid_x[2] )
    f2 = eval_LF1d_p( 2, L, grid_x, grid_x[2] )
    println("\nShould be close to zero: ", f1)
    println("\nShould be close to one: ", f2*sqrt(Δ))

    T = gen_T_matrix(L, grid_x)
    print("\nShould be close to zero: ")
    println(sum(abs.(T'-inv(T)))/Nbasis^2)

    ex_coefs = Array{Float64}(Nbasis)
    for i = 1:Nbasis
        ex_coefs[i] = myfunc(L, grid_x[i])
    end

    x = 2.0
    println("\nShould be close to each other:")
    println( myfunc(L, x) )
    println( eval_from_ex_coefs(ex_coefs*sqrt(Δ), L, grid_x, x) )

    FBR_ex_coefs = inv(T)*ex_coefs*sqrt(Δ)
    println( FBR_eval_from_ex_coefs(FBR_ex_coefs, L, x) )

    # FBR_ex_coefs = ifft(ex_coefs)
    # println( FBR_eval_from_ex_coefs(FBR_ex_coefs, L, x) )

end

main(10.0, 30)
