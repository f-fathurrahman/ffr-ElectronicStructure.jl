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


LATEX_START = """
\\documentclass[fleqn,a3paper,9pt]{article}
\\usepackage[a3paper]{geometry}
\\geometry{verbose,tmargin=2cm,bmargin=2cm,lmargin=2cm,rmargin=2cm}
\\usepackage{amsmath}

\\begin{document}

\\pagestyle{empty}

{\\footnotesize
"""

LATEX_END = """
}
\\end{document}
"""

HTML_START = """
<!DOCTYPE html>
<html>
<head>
<title>Laplacian Matrix</title>

<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
<meta http-equiv="X-UA-Compatible" content="IE=edge" />
<meta name="viewport" content="width=device-width, initial-scale=1">

<script type="text/x-mathjax-config">
  MathJax.Hub.Config({
    tex2jax: {inlineMath: [["$","$"],]}
  });
</script>
<script type="text/javascript"
  src="/home/efefer/mysoftwares/MathJax-master/MathJax.js?config=TeX-AMS_HTML-full">
</script>

</head>
<body>
"""

HTML_END = """
</body>
</html>
"""


# Create Laplacian matrix and its components
def create_nabla2( Nx, Ny, Nz ):

    DDX = createMatrix( "x", Nx )
    DDY = createMatrix( "y", Ny )
    DDZ = createMatrix( "z", Nz )

    # Identity matrices
    Ix = eye(Nx)
    Iy = eye(Ny)
    Iz = eye(Nz)

    tmp = do_kron(DDX,Iy)
    nabla2x = do_kron(tmp,Iz)

    tmp = do_kron(Ix,DDY)
    nabla2y = do_kron(tmp,Iz)

    tmp = do_kron(Ix,Iy)
    nabla2z = do_kron(tmp,DDZ)

    return nabla2x, nabla2y, nabla2z





def output_latex( Nx, Ny, Nz, nabla2x, nabla2y, nabla2z ):

    filname = "Laplacian_" + str(Nx) + "_" + str(Ny) + "_" + str(Nz) + ".tex"

    f = open(filname, "w")
    f.writelines(LATEX_START)

    f.write("\n")
    f.write("\\newpage\n")
    f.write("\\begin{equation*}\n")
    f.write(latex(nabla2x))
    f.write("\\end{equation*}\n")
    f.write("\n")

    f.write("\n")
    f.write("\\newpage\n")
    f.write("\\begin{equation*}\n")
    f.write(latex(nabla2y))
    f.write("\\end{equation*}\n")
    f.write("\n")

    f.write("\n")
    f.write("\\newpage\n")
    f.write("\\begin{equation*}\n")
    f.write(latex(nabla2z))
    f.write("\\end{equation*}\n")
    f.write("\n")

    f.write("\n")
    f.write("\\newpage\n")
    f.write("\\begin{equation*}\n")
    f.write(latex(nabla2x + nabla2y + nabla2z))
    f.write("\\end{equation*}\n")
    f.write("\n")

    f.writelines(LATEX_END)

    f.close()

def output_html( Nx, Ny, Nz, nabla2x, nabla2y, nabla2z ):

    filname = "Laplacian_" + str(Nx) + "_" + str(Ny) + "_" + str(Nz) + ".html"

    f = open(filname, "w")
    f.writelines(HTML_START)

    f.write("\n")
    f.write("\\begin{equation*}\n")
    f.write(latex(nabla2x))
    f.write("\\end{equation*}\n")
    f.write("\n")

    f.write("\n")
    f.write("\\begin{equation*}\n")
    f.write(latex(nabla2y))
    f.write("\\end{equation*}\n")
    f.write("\n")

    f.write("\n")
    f.write("\\begin{equation*}\n")
    f.write(latex(nabla2z))
    f.write("\\end{equation*}\n")
    f.write("\n")

    f.write("\n")
    f.write("\\begin{equation*}\n")
    f.write(latex(nabla2x + nabla2y + nabla2z))
    f.write("\\end{equation*}\n")
    f.write("\n")

    f.writelines(HTML_END)

    f.close()



import sys

Nargs = len(sys.argv)

if ( Nargs < 4 ):
    print("Error: this script needs 3 integer numbers as arguments")
    exit()

Nx = int( sys.argv[1] )
Ny = int( sys.argv[2] )
Nz = int( sys.argv[3] )

nablax, nablay, nablaz = create_nabla2( Nz, Ny, Nx )
output_html( Nx, Ny, Nz, nablax, nablay, nablaz )
