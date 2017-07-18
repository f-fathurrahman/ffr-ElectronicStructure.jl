using PyPlot
const plt = PyPlot

function plot_band_structure( evals, kpath; filename="BANDS.pdf" )

  Nkpts = size(kpath)[2]
  Nstates = size(evals)[1]

  Xcoords = zeros(Float64,Nkpts)
  Xcoords[1] = 0.0
  for ik = 2:Nkpts
    dk = kpath[:,ik] - kpath[:,ik-1]
    Xcoords[ik] = Xcoords[ik-1] + norm( dk )
  end

  plt.clf()
  for is = 1:Nstates
    y = evals[is,:]
    plt.plot( Xcoords[:], y, marker="o")
  end
  plt.grid()
  plt.savefig(filename)
end
