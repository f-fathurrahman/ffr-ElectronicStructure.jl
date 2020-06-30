# !INPUT/OUTPUT PARAMETERS:
#   z : atomic number (in,real)
# !DESCRIPTION:
#   Computes an approximate nuclear charge radius from the atomic number $Z$.
#   The nuclear mass number, $A$, is estimated using
#   $$ A=4.467\times 10^{-3}Z^2+2.163 Z-1.168, $$
#   [D. Andrae in {\it Relativistic Electronic Structure Theory - Fundamentals}
#   {\bf 11}, 203 (2002)], and the nuclear charge radius can be determined from
#   $$ r=\left(r_0+\frac{r_1}{A^{2/3}}+\frac{r_2}{A^{4/3}}\right)A^{1/3}, $$
#   where $r_0=0.9071$, $r_1=1.105$ and $r_2=-0.548$ [I. Angeli, {\it Atomic
#   Data and Nuclear Data Tables} {\bf 87}, 185 (2004)].

function radnucl( z::Float64 )
    
    # coefficients for computing mass number
    c2 = 4.467e-3
    c1 = 2.163
    c0 = -1.168
    
    # coefficients for computing charge radius (fm)
    r0 = 0.9071
    r1 = 1.105
    r2 = -0.548
    # Bohr radius in SI units (CODATA 2018)
    br_si = 0.529177210903e-10

    za = abs(z)
    a = 1.0
    # approximate nuclear mass number
    if za <= 1.0
        a = 1.0
    else
        a = c2*za^2 + c1*za + c0
    end
    # approximate nuclear charge radius
    a13 = a^(1.0/3.0)
    a23 = a13^2
    a43 = a13*a
    return (r0 + r1/a23 + r2/a43)*a13 * 1.0e-15/br_si

end

