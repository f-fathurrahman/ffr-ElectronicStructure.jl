"""
PyQuante in Julia
Experimenting with writing quantum chemistry in Julia

"""


function test()
  test_utils()
  test_pgbf()
  test_cgbf()
  test_overlap()
  test_kinetic()
  test_a_terms()
  test_gamma()
  test_na()
  test_fgamma()
  test_one()
  test_na2()
  test_two_terms()
  test_coul1()
  test_vrr()
  test_hrr()
  test_geo_basis()
  test_h2()
  test_lih()
  test_h2o()
end

test()
