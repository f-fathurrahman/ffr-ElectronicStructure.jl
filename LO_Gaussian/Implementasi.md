## Beberapa shortcut

```julia
const Tuple3F64 = Tuple{Float64,Float64,Float64}
const Tuple3I64 = Tuple{Int64,Int64,Int64}
```

## Fungsi Gaussian primitif

Tipe data: `PGBF`

```julia
struct PGBF
    center::Tuple3F64
    power::Tuple3I64
    expn::Float64
    NORM::Float64
end
```

`NORM` (semua dalam huruf kapital) untuk membedakan dari fungsi `norm`.
`expn` untuk membedakan dari fungsi matematik `exp`.
