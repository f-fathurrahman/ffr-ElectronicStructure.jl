factorial2(n) = prod(n:-2:1) # double factorial !!
dist2(dx,dy,dz) = dx*dx+dy*dy+dz*dz # Is there something in the standard library that does this?

function pairs(n::Int64,which="diag")
  function _it()
    for j in 1:n
      start = j+1
      if which=="diag"
        start = j
      elseif which=="rect"
        start = 1
      end
      for i in start:n
        produce((j,i))
      end
    end
  end
  Task(_it)
end

triangle(i::Int64) = div(i*(i+1),2)
triangle(i::Int64,j::Int64) = i<j ? triangle(j-1)+i : triangle(i-1)+j

function iiterator(n::Int64)
  function _it()
    for (i,j) in pairs(n)
      ij = triangle(i-1,j)
      for (k,l) in pairs(n)
        kl = triangle(k-1,l)
        if ij <= kl
          produce((i,j,k,l))
        end
      end
    end
  end
  Task(_it)
end

iindex(i::Int64,j::Int64,k::Int64,l::Int64) = triangle(triangle(i,j),triangle(k,l))

trace2(A,B) = sum(A.*B)

function test_utils()
  @assert factorial2(6)==48
  @assert collect(pairs(3)) == Any[(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]
  @assert collect(pairs(3,"subdiag")) == Any[(1,2),(1,3),(2,3)]
  @assert collect(pairs(2,"rect")) == Any[(1,1),(1,2),(2,1),(2,2)]
  @assert iindex(1,1,1,1) == 1
  @assert iindex(1,1,1,2) == iindex(1,1,2,1) == iindex(1,2,1,1) == iindex(2,1,1,1) == 2
  @assert iindex(1,1,2,2) == iindex(2,2,1,1) == 4
  @printf "test_utils is passed\n"
end
