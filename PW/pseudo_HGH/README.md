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

## Prototype of `op_V_ps_NL`

without k-dependence

Using notation in Kohanoff's instead of beta ?

f_NL includes structure factor.

```julia
for ig = 1:Ngwx
  F_NL[ig,slm,ist] = f_NL[ig,slm,] * psi[ig,ist]
end
```


slm: combined indices of species, lm, ...

probably also needs index of iproj


```julia
function op_V_ps_NL( pw, f_NL, F_NL, psi )

  Ngwx = size(psi)[1]
  Nstates = size(psi)[2]

  out = zeros( Complex128, Ngwx, Nstates )

  for ig = 1:Ngwx
    out[ig] =     
  end

end
```

