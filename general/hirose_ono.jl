"""
An attempt to understand Hirose-Ono paper entitled:

  Direct minimization to generate electronic states with proper
  occupation numbers
  Phys. Rev. B. 64 085105
"""

function test_main()

  #srand(1234)
  Nstates = 5
  Npoints = 10

  ϕ = rand(Npoints,Nstates) + im*rand(Npoints,Nstates)

  for ist = 2:Nstates
    d = dot( ϕ[:,ist], ϕ[:,1] )
    @printf("States: %4d d = [%18.10f,%18.10f]\n", ist, real(d), imag(d))
  end

  U = ϕ'*ϕ

  Um = inv(sqrtm(U))  # Um = U^{-1/2}

  ϕ_ortho = ϕ * Um

  # Check orthogonalization wrt 1st states
  println("")
  for ist = 2:Nstates
    d = dot( ϕ_ortho[:,ist], ϕ_ortho[:,1] )
    @printf("States: %4d d = [%18.10f,%18.10f]\n", ist, real(d), imag(d))
  end

  #c = copy(Um)
  c = sqrtm(U)   # <----- this is the c_{l,i}

  phi = zeros(Complex128,Npoints,Nstates)
  for i = 1:Nstates
    for l = 1:Nstates
      phi[:,i] = phi[:,i] + c[l,i] * ϕ_ortho[:,l]
    end
  end

  # Check orthogonalization wrt 1st states
  println("")
  for ist = 2:Nstates
    d = dot( phi[:,ist], phi[:,1] )
    @printf("States: %4d d = [%18.10f,%18.10f]\n", ist, real(d), imag(d))
  end

  T = zeros(Complex128,Nstates,Nstates)
  for k = 1:Nstates
    for l = 1:Nstates
      for i = 1:Nstates
        T[k,l] = T[k,l] + c[k,i]*c[l,i]
      end
    end
  end

  λ = eigvals(T)
  Focc = real((2 - λ).*λ)  # the occupation number, should not exceed 1.0

  println("\nOccupation numbers")
  for ist = 1:Nstates
    @printf("%4d %18.10f\n", ist, Focc[ist])
  end

end

test_main()
