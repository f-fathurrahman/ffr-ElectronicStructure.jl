The HGH pseudopotentials used are given in ABINIT format.

NL pseudopotential evaluation based on PWSCF.

Evaluating spherical harmonics

```julia
# Get lmaxkb
lmaxkb = 3  # or use PsNL_lmax ?

# get active G-vectors for wavefunction
idx_gw2r = pw.gvecw.idx_gw2r
Gk = pw.gvec.G[:,idx_gw2r] + kpt[:]

ylm = ylmr2( lmaxkb, Gk )
```


Evaluating beta functions (`init_us_2`)

```julia
for isp = 1:Nspecies

  for i = 1:PsPot[isp].snprj

     l = PsPot[isp].lll[i]
     iproj = PsPot[isp].ipr[i]
     Vprj = eval_HGH_proj_G( psp[isp], l, iproj, G, Ω )

  end

end
```
