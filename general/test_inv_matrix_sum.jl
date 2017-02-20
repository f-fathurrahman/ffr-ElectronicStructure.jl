function test_main()

  srand(1234)

  A = rand(3,3)
  A = 0.5*(A + A')
  invA = inv(A)

  println("\nMatrix A:")
  println(A)
  println("\nInverse of A")
  println(invA)

  B = rand(3,3)
  B = 0.5*(B + B')
  invB = inv(B)

  println("\nMatrix B")
  println(B)
  println("\nInverse of B")
  println(invB)

  AB = A + B
  invAB = inv(A + B)

  println("\nMatrix A + B")
  println(AB)
  println("\nInverse of A + B")
  println(invAB)

  test_inv_Miller( A, B )

  println(rank(B))

end

function test_inv_Miller( G, H )

  g = trace( H*inv(G) )
  println("g = ", g)

  invGH_mill = inv(G) - 1/(1+g) * inv(G)*H*inv(G)

  println("\nTest Miller")
  println(invGH_mill - inv(G+H))
end


test_main()
