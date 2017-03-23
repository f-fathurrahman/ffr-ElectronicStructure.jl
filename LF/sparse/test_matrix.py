from sympy import *

Nx = 3
Ny = 3
Nz = 3

# Symbols for matrix elements
x_11, x_12, x_13 = symbols("x_11 x_12 x_13")
x_21, x_22, x_23 = symbols("x_21 x_22 x_23")
x_31, x_32, x_33 = symbols("x_31 x_32 x_33")

y_11, y_12, y_13 = symbols("y_11 y_12 y_13")
y_21, y_22, y_23 = symbols("y_21 y_22 y_23")
y_31, y_32, y_33 = symbols("y_31 y_32 y_33")

z_11, z_12, z_13 = symbols("z_11 z_12 z_13")
z_21, z_22, z_23 = symbols("z_21 z_22 z_23")
z_31, z_32, z_33 = symbols("z_31 z_32 z_33")

# Identity matrices
Ix = eye(Nx)
Iy = eye(Ny)
Iz = eye(Nz)


DDX = Matrix([ [x_11, x_12, x_13], [x_21, x_22, x_23], [x_31, x_32, x_33] ])
DDY = Matrix([ [y_11, y_12, y_13], [y_21, y_22, y_23], [y_31, y_32, y_33] ])
DDZ = Matrix([ [z_11, z_12, z_13], [z_21, z_22, z_23], [z_31, z_32, z_33] ])

# kronecker product between DDX and Iy
#tmpxy = zeros( Nx+Ny, Nx+Ny )
#tmpxy[0:3,0:3] = DDX[0,0]*Iy
#tmpxy[0:3,3:Nx+Ny] = DDX[0,1]*Iy
#tmpxy[3:Nx+Ny,0:3] = DDX[1,0]*Iy
#tmpxy[3:Nx+Ny,3:Nx+Ny] = DDX[1,1]*Iy

# A and B should be square matrices
def do_kron( A, B ):
    nA = A.shape[0]
    nB = B.shape[0]
    Z = zeros(nA*nB,nA*nB)
    #print("Shape of Z:", Z.shape)
    # loop over cols
    for j in range(0,nA):
        jj = j*nB
        #print("")
        # loop over rows
        for i in range(0,nA):
            ii = i*nB
            #print("Range rows: %d %d" % (ii, ii+nB) )
            #print("Range cols: %d %d" % (jj, jj+nB) )
            Z[ii:ii+nB,jj:jj+nB] = A[i,j]*B
    #
    return Z



tmp = do_kron(DDX,Iy)
nablax = do_kron(tmp,Iz)
print(latex(nablax))

tmp = do_kron(Ix,DDY)
nablay = do_kron(tmp,Iz)
print(latex(nablay))

tmp = do_kron(Ix,Iy)
nablaz = do_kron(tmp,DDZ)
print(latex(nablaz))
