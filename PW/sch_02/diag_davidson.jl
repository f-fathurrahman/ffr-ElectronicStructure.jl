function diag_davidson( pw::PWGrid, Vpot, X0;
                        tol=1e-5, tol_avg=1e-7, NiterMax=200, verbose=false )

  # get size info
  Ncols = size(X0)[2]
  Ngwx  = size(X0)[1]

  if Ncols <= 0
   @printf("diag_davidson requires at least one initial wave function!\n");
   return
  end
  # orthonormalize the initial wave functions.
  X = ortho_gram_schmidt(X0)  # normalize (again)?

  evals = zeros(Float64,Ncols)
  R = zeros(Complex128,Ngwx,Ncols)
  Hred = zeros(Complex128, 2*Ncols, 2*Ncols)
  Sred = zeros(Complex128, 2*Ncols, 2*Ncols)
  resnrm = zeros(Float64,Ncols)

  HX = apply_H( pw, Vpot, X )

  for i = 1:Ncols
    evals[i] = real( dot(X[i], HX[i]) )
  end

  for i = 1:Ncols
    R[:,i] = evals[i]*X[:,i] - HX[:,i]
    resnrm[i] = real( dot( R[:,i], R[:,i] ) )
  end

  for iter = 1:10

    R = Kprec(pw, R)

    HR = apply_H( pw, Vpot, R )

    # FIXME: Pull this outside the loop?
    if iter == 1
      Hred[1:Ncols,1:Ncols] = X' * HX
    else
      Hred[1:Ncols,1:Ncols] = diagm(evals)
    end
    Hred[1:Ncols,Ncols+1:2*Ncols] = X' * HR
    Hred[Ncols+1:2*Ncols,Ncols+1:2*Ncols] = R' * HR
    Hred[Ncols+1:2*Ncols,1:Ncols] = Hred[1:Ncols,Ncols+1:2*Ncols]'

    Sred[1:Ncols,1:Ncols] = X' * X
    Sred[1:Ncols,Ncols+1:2*Ncols] = X' * R
    Sred[Ncols+1:2*Ncols,Ncols+1:2*Ncols] = R' * R
    Sred[Ncols+1:2*Ncols,1:Ncols] = Sred[1:Ncols,Ncols+1:2*Ncols]'

    λ_red, X_red = eig(Hred,Sred)

    evals[:] = real(λ_red[1:Ncols])

    X = X * X_red[1:Ncols,1:Ncols] + R*X_red[Ncols+1:2*Ncols,1:Ncols]
    HX = HX * X_red[1:Ncols,1:Ncols] + HR*X_red[Ncols+1:2*Ncols,1:Ncols]

    X = ortho_gram_schmidt(X)

    @printf("\n")
    for i = 1:Ncols
      R[:,i] = evals[i]*X[:,i] - HX[:,i]
      @printf("%4d %18.10f %18.10f\n", i, evals[i], norm(R[:,i]))
    end

  end

  exit()

end
