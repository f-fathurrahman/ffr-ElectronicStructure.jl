function op_H( pw, Vpot, psi )
    return op_K( pw, psi ) + op_Vpot( pw, Vpot, psi )
end
