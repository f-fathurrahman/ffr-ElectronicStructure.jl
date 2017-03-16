const XC_LDA_X = 1   # Exchange
const XC_LDA_C_VWN = 7  # Vosko, Wilk, & Nusair (5)
const XC_GGA_X_PBE = 101
const XC_GGA_C_PBE = 130

function test_main()

  ρ = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
  σ = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
  N = size(ρ)[1]

  ϵ_x = zeros(N)
  ϵ_c = zeros(N)

  ccall( (:ffrlibxc, "ffrlibxc"), Void,
         (Int, Int, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
         XC_LDA_X, N, ρ, σ, ϵ_x )

  ccall( (:ffrlibxc, "ffrlibxc"), Void,
         (Int, Int, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          XC_GGA_C_PBE, N, ρ, σ, ϵ_c )

  for i = 1:N
    @printf("%d %f %f %f\n", i, ρ[i], σ[i], ϵ_x[i] + ϵ_c[i])
  end

end

test_main()
