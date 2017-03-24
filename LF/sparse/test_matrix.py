from sympy import *

def createMatrix( symb, N ):
    rows = []
    for i in range(1,N+1):
        str1 = symb + "_" + str(i) + "1:" + str(i) + str(N+1)
        #print str1
        rows.append( list( symbols(str1) ) )
    #
    return Matrix(rows)


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


STR_START = """
\\documentclass[fleqn,a3paper,9pt]{article}
\\usepackage[a3paper]{geometry}
\\geometry{verbose,tmargin=2cm,bmargin=2cm,lmargin=2cm,rmargin=2cm}
\\usepackage{amsmath}

\\begin{document}

\\pagestyle{empty}

{\\footnotesize
"""

STR_END = """
}
\\end{document}
"""

def do_latex( Nx, Ny, Nz ):

    DDX = createMatrix( "x", Nx )
    DDY = createMatrix( "y", Ny )
    DDZ = createMatrix( "z", Nz )

    # Identity matrices
    Ix = eye(Nx)
    Iy = eye(Ny)
    Iz = eye(Nz)

    tmp = do_kron(DDX,Iy)
    nablax = do_kron(tmp,Iz)

    print(STR_START)

    print("")
    print("\\newpage")
    print("\\begin{equation*}")
    print(latex(nablax))
    print("\\end{equation*}")
    print("")

    tmp = do_kron(Ix,DDY)
    nablay = do_kron(tmp,Iz)
    print("")
    print("\\newpage")
    print("\\begin{equation*}")
    print(latex(nablay))
    print("\\end{equation*}")
    print("")

    tmp = do_kron(Ix,Iy)
    nablaz = do_kron(tmp,DDZ)
    print("")
    print("\\newpage")
    print("\\begin{equation*}")
    print(latex(nablaz))
    print("\\end{equation*}")
    print("")

    print("")
    print("\\newpage")
    print("\\begin{equation*}")
    print(latex(nablax + nablay + nablaz))
    print("\\end{equation*}")
    print("")

    print(STR_END)


do_latex( 3, 3, 3 )
