function read_kpts( filname )
  file = open(filname)
  str = readline(file)
  Nkpts = parse( Int, str )
  kpts = zeros( Float64, 3,Nkpts )
  for ik = 1:Nkpts
    str = split(readline(file))
    kpts[1,ik] = parse( Float64, str[1] )
    kpts[2,ik] = parse( Float64, str[2] )
    kpts[3,ik] = parse( Float64, str[3] )
  end
  close(file)
  return kpts
end
