const ⊗ = kron

function test_main()

  srand(1234)
  A = rand(3,3); A = 0.5*(A+A')
  B = rand(3,3); B = 0.5*(B+B')
  C = rand(3,3); C = 0.5*(C+C')

  term1 = A ⊗ eye(3) ⊗ eye(3)
  term2 = eye(3) ⊗ B ⊗ eye(3)
  term3 = eye(3) ⊗ eye(3) ⊗ C

  Lapl = term1 + term2 + term3

  invLapl = inv(Lapl)

  invterm1 = inv(term1)
  invterm2 = inv(term2)
  invterm3 = inv(term3)

  # these should be zeros
  diff1 = invterm1 - inv(A) ⊗ eye(3) ⊗ eye(3)
  diff2 = invterm2 - eye(3) ⊗ inv(B) ⊗ eye(3)
  diff3 = invterm3 - eye(3) ⊗ eye(3) ⊗ inv(C)

  # TODO: Find relationship between invterm1 and invLapl

end
