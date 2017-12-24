function op_H( pw, Potentials, psi )
    V_loc = Potentials.Ionic + Potentials.Hartree + Potentials.XC
    return op_K( pw, psi ) + op_V_loc( pw, V_loc, psi )
end
