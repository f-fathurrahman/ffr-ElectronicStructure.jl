The HGH pseudopotentials used are given in ABINIT format.

NL pseudopotential evaluation based on PWSCF.

Evaluating spherical harmonics

```julia
# Get lmaxkb
lmaxkb = 3  # or use PsNL_lmax ?

```


Evaluating beta functions (`init_us_2`)

```julia
for isp = 1:Nspecies

  for ibeta = 1:PsPot[isp].nbeta

     l = PsPot.lll[ibeta]
     iproj = PsPot[isp].ipr[ibeta]
     Vprj = eval_HGH_proj_G( psp[isp], l, iproj, G, Ω )

  end

end
```
