struct LF3dGrid
  LFx::LF1dGrid
  LFy::LF1dGrid
  LFz::LF1dGrid
  #
  Nx::Int64
  Ny::Int64
  Nz::Int64
  # XXX probably these variables are not important
  Lx::Float64
  Ly::Float64
  Lz::Float64
  #
  lingrid::Array{Float64,2}
  xyz2lin::Array{Int64,3}
  lin2xyz::Array{Int64,2}
end


function init_LF3d_p( NN::Array{Int64,1}, AA::Array{Float64,1}, BB::Array{Float64,1};
                      verbose=false )
  #
  Nx = NN[1]
  Ny = NN[2]
  Nz = NN[3]
  #
  Lx = BB[1] - AA[1]
  Ly = BB[2] - AA[2]
  Lz = BB[3] - AA[3]
  #
  LFx = init_LF1d_p( Nx, AA[1], BB[1], verbose )
  LFy = init_LF1d_p( Ny, AA[2], BB[2], verbose )
  LFz = init_LF1d_p( Nz, AA[3], BB[3], verbose )
  #
  lingrid = zeros(Float64, 3, Nx*Ny*Nz )
  xyz2lin = zeros(Int64, Nx, Ny, Nz )
  lin2xyz = zeros(Int64, 3, Nx*Nz*Ny )
  ip = 0
  for k = 1 : Nz
    for j = 1 : Ny
      for i = 1 : Nx
        ip = ip + 1
        #
        lingrid[1,ip] = LFx.grid[i]
        lingrid[2,ip] = LFy.grid[j]
        lingrid[3,ip] = LFz.grid[k]
        #
        xyz2lin[i,j,k] = ip
        lin2xyz[1:3,ip] = [i,j,k]
      end
    end
  end
  if verbose
    @printf("Initializing an instance of LF3dGrid: periodic in xyz directions.\n")
  end
  LF = LF3dGrid( LFx,LFy,LFz, Nx,Ny,Nz, Lx,Ly,Lz, lingrid, xyz2lin, lin2xyz )
  return LF
end


function init_LF3d_c( NN::Array{Int64,1}, AA::Array{Float64,1}, BB::Array{Float64,1};
                      verbose=false )
  #
  Nx = NN[1]
  Ny = NN[2]
  Nz = NN[3]
  #
  Lx = BB[1] - AA[1]
  Ly = BB[2] - AA[2]
  Lz = BB[3] - AA[3]
  #
  LFx = init_LF1d_c( Nx, AA[1], BB[1] )
  LFy = init_LF1d_c( Ny, AA[2], BB[2] )
  LFz = init_LF1d_c( Nz, AA[3], BB[3] )
  #
  lingrid = zeros(Float64, 3, Nx*Ny*Nz )
  xyz2lin = zeros(Int64, Nx, Ny, Nz )
  lin2xyz = zeros(Int64, 3, Nx*Nz*Ny )
  ip = 0
  for k = 1 : Nz
    for j = 1 : Ny
      for i = 1 : Nx
        ip = ip + 1
        #
        lingrid[1,ip] = LFx.grid[i]
        lingrid[2,ip] = LFy.grid[j]
        lingrid[3,ip] = LFz.grid[k]
        #
        xyz2lin[i,j,k] = ip
        lin2xyz[1:3,ip] = [i,j,k]
      end
    end
  end
  if verbose
    @printf("Initializing an instance of LF3dGrid: cluster BC in xyz directions.\n")
  end
  LF = LF3dGrid( LFx,LFy,LFz, Nx,Ny,Nz, Lx,Ly,Lz, lingrid, xyz2lin, lin2xyz )
  return LF
end


function init_LF3d_sinc( NN::Array{Int64,1}, hh::Array{Float64,1}; verbose=false )
  #
  Nx = NN[1]
  Ny = NN[2]
  Nz = NN[3]
  #
  LFx = init_LF1d_sinc( Nx, hh[1] )
  LFy = init_LF1d_sinc( Ny, hh[2] )
  LFz = init_LF1d_sinc( Nz, hh[3] )
  #
  Lx = LFx.B - LFy.A
  Ly = LFy.B - LFy.A
  Lz = LFz.B - LFz.A
  #
  lingrid = zeros(Float64, 3, Nx*Ny*Nz )
  xyz2lin = zeros(Int64, Nx, Ny, Nz )
  lin2xyz = zeros(Int64, 3, Nx*Nz*Ny )
  ip = 0
  for k = 1 : Nz
    for j = 1 : Ny
      for i = 1 : Nx
        ip = ip + 1
        #
        lingrid[1,ip] = LFx.grid[i]
        lingrid[2,ip] = LFy.grid[j]
        lingrid[3,ip] = LFz.grid[k]
        #
        xyz2lin[i,j,k] = ip
        lin2xyz[1:3,ip] = [i,j,k]
      end
    end
  end
  if verbose
    @printf("Initializing an instance of LF3dGrid: Lagrange-sinc in xyz directions.\n")
  end
  LF = LF3dGrid( LFx,LFy,LFz, Nx,Ny,Nz, Lx,Ly,Lz, lingrid, xyz2lin, lin2xyz )
  return LF
end
