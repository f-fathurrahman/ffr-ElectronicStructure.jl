using LinearAlgebra

struct LatticeVars
    # lattice vectors stored column-wise
    avec::Array{Float64,2}
    # magnitude of random displacements added to lattice vectors
    rndavec::Float64
    # inverse of lattice vector matrix
    ainv::Array{Float64,2}
    # reciprocal lattice vectors
    bvec::Array{Float64,2}
    # inverse of reciprocal lattice vector matrix
    binv::Array{Float64,2}
    # unit cell volume
    omega::Float64
    # Brillouin zone volume
    omegabz::Float64
    # any vector with length less than epslat is considered zero
    epslat::Float64
end

#=
case('avec')
  read(50,*,err=20) avec(:,1)
  read(50,*,err=20) avec(:,2)
  read(50,*,err=20) avec(:,3)

avec is stored like lattice vectors in Quantum Espresso
=#

function LatticeVars( LatVecs )
    @assert size(LatVecs,1) == 3
    @assert size(LatVecs,2) == 3
    
    avec = LatVecs[:,:]
    ainv = inv(avec)
    
    bvec = 2*pi*inv(avec')
    binv = inv(bvec)
    
    omega = det(LatVecs)
    omegabz = (2*pi)^3/omega

    epslat = 1e-6
    rndavec = 0.0

    return LatticeVars( avec, rndavec, ainv, bvec, binv, omega, omegabz, epslat )
end