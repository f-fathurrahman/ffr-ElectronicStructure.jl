# double factorial
factorial2(n) = prod(n:-2:1)

dist2(dx,dy,dz) = dx*dx + dy*dy + dz*dz

triangle(i::Int64) = div(i*(i+1),2)
triangle(i::Int64,j::Int64) = i < j ? triangle(j-1)+i : triangle(i-1)+j

iindex(i::Int64,j::Int64,k::Int64,l::Int64) = triangle(triangle(i,j),triangle(k,l))

trace2(A,B) = sum(A.*B)
