# ffr, 19 April 2016

module m_LF1d

# for the moment, I decided to not include L, because it can
# be calculated from A and B
type LF1dGrid
  N::Int64
  A::Float64
  B::Float64
  h::Float64
  grid::Array{Float64,1}
  D1jl::Array{Float64,2}
  D2jl::Array{Float64,2}
end

# Periodic LF
function init_LF1d_p( N::Int64, A::Float64, B::Float64, verbose=false )
  # Check argument
  if N % 2 == 0
    error("N must be an odd number\n")
  end
  L = (B - A)
  h = L/N
  grid = zeros(Float64,N)
  for i = 1:N
    grid[i] = A + 0.5*(B-A)*(2*i-1)/N
  end
  D1jl = zeros(Float64,N,N)
  D2jl = zeros(Float64,N,N)
  #
  # Diagonal elements
  #
  Nprimed = (N-1)/2
  for j = 1:N
    # the diagonal elements of D1jl are already zero
    D2jl[j,j] = -(2.*pi/L)^2 * Nprimed * (Nprimed+1)/3.
  end
  #
  # Off diagonal elements
  #
  for j = 1 : N
    for l = j+1 : N
      #
      nn = j - l
      #
      tt1 = pi/L * (-1.0)^nn
      tt2 = sin(pi*nn/N)
      #
      tt3 = (2.0*pi/L)^2 * (-1.0)^nn * cos(pi*nn/N)
      tt4 = 2.0*sin(pi*nn/N)^2
      #
      D1jl[j,l] =  tt1/tt2
      D1jl[l,j] = -tt1/tt2
      #
      D2jl[j,l] = -tt3/tt4
      D2jl[l,j] = -tt3/tt4
    end
  end
  #
  LF = LF1dGrid( N, A, B, h, grid, D1jl, D2jl )
  if verbose
    @printf("Allocated: periodic 1d LBF grid: N, h = %d , %f\n", N, h)
  end
  #
  return LF
end # function init_LF1d_p


# cluster Lagrange function
function init_LF1d_c( N::Int64, A::Float64, B::Float64, verbose=false )
  L = B - A
  h = L/(N + 1.)
  #
  grid = zeros(Float64,N)
  for i = 1:N
    grid[i] = A + i*(B-A)/(N+1)
  end
  #
  D1jl = zeros(Float64,N,N)  # XXX this is not yet used
  D2jl = zeros(Float64,N,N)
  #
  # Diagonal part
  #
  pre = -pi^2/(2.*L^2)
  for i = 1:N
    t1 = ( 2.0*(N+1)^2 + 1 )/3.0
    t2 = sin( i*pi/(N+1) )^2
    D2jl[i,i] = pre*( t1 - 1.0/t2 )
  end
  #
  # Off-diagonal
  #
  for l = 1 : N
    for j = l+1 : N
      nnm = l - j
      nnp = l + j
      pre = -pi^2 / (2*L^2) * (-1.0)^nnm
      t1 = sin( pi*nnm/2./(N+1) )^2
      t2 = sin( pi*nnp/2./(N+1) )^2
      #
      D2jl[l,j] = pre*( 1.0/t1 - 1.0/t2 )
      D2jl[j,l] = pre*( 1.0/t1 - 1.0/t2 )  # XXX is it faster to just recalculate?
    end
  end
  LF = LF1dGrid( N, A, B, h, grid, D1jl, D2jl )
  if verbose
    @printf("Allocated: cluster 1d LBF grid: N, h = %5d , %10.5f\n", N, h)
  end
  return LF
end


# Lagrange-sinc function
function init_LF1d_sinc( N::Int64, h::Float64, verbose=false )
  # Choice for A and B
  A = -(N-1)*h/2.0
  B =  (N-1)*h/2.0
  #
  grid = zeros(Float64,N)
  for i = 1:N
    grid[i] = A + (i-1)*h
  end
  #
  D1jl = zeros(Float64,N,N)  # XXX this is not yet used
  D2jl = zeros(Float64,N,N)
  #
  # Diagonal part
  #
  for i = 1:N
    D2jl[i,i] = -pi^2 / 3.0 / h^2
  end
  #
  # Off-diagonal
  #
  for j = 1 : N
    for i = j+1 : N
      D2jl[i,j] = -2.0*(-1.0)^(i-j)/( grid[i] - grid[j] )^2
      D2jl[j,i] = D2jl[i,j]  # XXX is it faster to just recalculate?
    end
  end
  LF = LF1dGrid( N, A, B, h, grid, D1jl, D2jl )
  if verbose
    @printf("Allocated: 1d Lagrange-sinc grid: N, h = %5d , %10.5f\n", N,  h)
  end
  return LF
end


export LF1dGrid
export init_LF1d_p
export init_LF1d_c
export init_LF1d_sinc

end # module m_LF1d

