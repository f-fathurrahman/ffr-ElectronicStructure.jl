include("spharm.jl")

function test_main()

  Ng = 5
  lmax = 2

  Gvec = zeros(3,Ng)

  for ig = 1:Ng
    Gvec[:,ig] = [1.3, 1.1, 2.1]*(ig-1)
  end

  ylm = zeros( Float64, Ng, (lmax+1)^2 )

  ilm = 0
  for l = 0:lmax
    @printf("l = %d\n", l)
    # m = 0
    ilm = ilm + 1
    ylm[:,ilm] = spharm( l, 0, Gvec )
    println(ilm)
    println(ilm, " ", l, " ", 0)
    # m > 0 and m < 0
    for m = 1:l
      ilm = ilm + 1
      ylm[:,ilm] = spharm( l, m, Gvec )*(-1.0)^m
      println(ilm, " ", l, " ", m)
      #
      ilm = ilm + 1
      ylm[:,ilm] = spharm( l, -m, Gvec )*(-1.0)^m
      println(ilm)
      println(ilm, " ", l, " ", -m)
    end
  end

  for lm = 1:(lmax+1)^2
    @printf("\n")
    @printf("lm = %5d\n", lm)
    for ig = 1:Ng
      @printf("%5d %18.10f\n", ig, ylm[ig,lm])
    end
  end

end

test_main()
