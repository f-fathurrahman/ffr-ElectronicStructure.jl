function ILU(v)
  n = size(v,1)
  e = ones(n,1)
  L = spdiags( [-e e], -1:0, n, n )
  U = L'
  return U\(L\v)
end

#------------- Build a diagonal matrix
function spdiags(B,d,m,n)
  #   spdiags(B,d,m,n)
  # creates a sparse matrix from its diagonals

  d = d[:]
  p = length(d)

  len = zeros(p+1,1)
  for k = 1:p
    len[k+1] = round( Int, len[k]+length(max(1,1-d[k]):min(m,n-d[k])) )
  end
  a = zeros(round(Int,len[p+1]),3)
  for k = 1:p
    # Append new d[k]-th diagonal to compact form
    i = max(1,1-d[k]):min(m,n-d[k])
    a[(round(Int,len[k])+1):round(Int,len[k+1]),:] = [i i+d[k] B[i+(m>=n)*d[k],k]]
  end

  A = sparse(round(Int,a[:,1]),round(Int,a[:,2]),a[:,3],m,n);

  return A
end
