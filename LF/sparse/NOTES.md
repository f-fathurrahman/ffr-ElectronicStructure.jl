## Sparse matrix in Julia

Type `SparseMatrixCSC` has 5 fields: `m`, `n`, `nzval`, `rowval`
and `colptr`.

Create a random sparse matrix:

```julia
julia> A = sprand(4,4,0.4)
4x4 sparse matrix with 7 Float64 entries:
  [2, 2]  =  0.316691
  [1, 3]  =  0.504048
  [2, 3]  =  0.985886
  [3, 3]  =  0.328337
  [4, 3]  =  0.467175
  [1, 4]  =  0.969968
  [2, 4]  =  0.0916004

julia> full(A)
4x4 Array{Float64,2}:
  0.0  0.0       0.504048  0.969968
  0.0  0.316691  0.985886  0.0916004
  0.0  0.0       0.328337  0.0      
  0.0  0.0       0.467175  0.0
```

`m` and `n` are the dimensions of the matrix (m: number of rows, n: number
of columns).


colptr starts from 1. Then, as we loop through rows in the first and count
number of non-zeros, we add this number to the first element of colptr.
When we found certain row with non-zero value, we add row index to
rowval.
As an example for the matrix A above, we don't have any non-zero element
for the first column, then the second element of colptr is the same, i.e 1.
So, after looping through elements in the first column we have:

```julia
colptr = [1, 1]
rowval = []
nzval = []
```

The for the second column, we have one non-zero element at row 2, so now
we have

```julia
colptr = [1, 1, 2]
rowval = [2]
nzval = [0.316691]
```

Similarly, after the 3rd column (we have 4 non-zero values)

```julia
colptr = [1, 1, 2, 6]
rowval = [2, 1, 2, 3, 4]
nzval = [0.316691, 0.504048 0.985886, 0.328337, 0.467175]
```

Finally after 4th column (we have 2 non-zero values)

```julia
colptr = [1, 1, 2, 6, 8]
rowval = [2, 1, 2, 3, 4, 1, 2]
nzval = [0.316691, 0.504048 0.985886, 0.328337, 0.467175, 0.969968, 0.0916004]
```
