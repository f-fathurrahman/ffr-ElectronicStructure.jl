using Printf

include("ylmr2.jl")

function test_main()

  Ng = 5
  lmax = 3

  Gvec = zeros(3,Ng)

  for ig = 1:Ng
    Gvec[:,ig] = [-2.1, 2.0, 1.0]
  end

  ylm = ylmr2( lmax, Gvec )

  for lm = 1:(lmax+1)^2
    @printf("\n")
    @printf("lm = %5d\n", lm)
    for ig = 1:Ng
      @printf("%5d %18.10f\n", ig, ylm[ig,lm])
    end
  end

end

test_main()
