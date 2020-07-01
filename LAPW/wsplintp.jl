function wsplintp!(n::Int64, x, w)
    ## arguments
    #integer, intent(in) :: n
    #real(8), intent(in) :: x(n)
    #real(8), intent(out) :: w(4,n)
    
    @assert n >= 4
    
    w[:,1] .= 0.0
    f = zeros(4)

    f[1] = 1.0
    w[1,2] = polynm(-1, 4, x, f, x[2])
    
    f[1] = 0.0
    f[2] = 1.0
    w[2,2] = polynm(-1, 4, x, f, x[2])
    f[2] = 0.0
    
    f[3] = 1.0
    w[3,2] = polynm(-1, 4, x, f, x[2])
    f[3] = 0.0
    
    f[4] = 1.0
    w[4,2] = polynm(-1, 4, x, f, x[2])

    for i in 3:n-1
        f[:] .= 0.0
        f[1] = 1.0
        @views t1 = polynm(-1, 4, x[i-2:end], f, x[i-1])
        @views t2 = polynm(-1, 4, x[i-2:end], f, x[i])
        w[1,i] = t2 - t1
        
        f[1] = 0.0
        f[2] = 1.0
        @views t1 = polynm(-1, 4, x[i-2:end], f, x[i-1])
        @views t2 = polynm(-1, 4, x[i-2:end], f, x[i])
        w[2,i] = t2 - t1
        
        f[2] = 0.0
        f[3] = 1.0
        @views t1 = polynm(-1, 4, x[i-2:end], f, x[i-1])
        @views t2 = polynm(-1, 4, x[i-2:end], f, x[i])
        w[3,i] = t2 - t1
        
        f[3] = 0.0
        f[4] = 1.0
        @views t1 = polynm(-1, 4, x[i-2:end], f, x[i-1])
        @views t2 = polynm(-1, 4, x[i-2:end], f, x[i])
        w[4,i] = t2 - t1
    end

    f[:] .= 0.0
    f[1] = 1.0
    @views t1 = polynm(-1, 4, x[n-3:end], f, x[n-1])
    @views t2 = polynm(-1, 4, x[n-3:end], f, x[n])
    w[1,n] = t2 - t1

    f[1] = 0.0
    f[2] = 1.0
    @views t1 = polynm(-1, 4, x[n-3:end], f, x[n-1])
    @views t2 = polynm(-1, 4, x[n-3:end], f, x[n])
    w[2,n] = t2 - t1
    
    f[2] = 0.0
    f[3] = 1.0
    @views t1 = polynm(-1, 4, x[n-3:end], f, x[n-1])
    @views t2 = polynm(-1, 4, x[n-3:end], f, x[n])
    w[3,n] = t2 - t1
    
    f[3] = 0.0
    f[4] = 1.0
    @views t1 = polynm(-1, 4, x[n-3:end], f, x[n-1])
    @views t2 = polynm(-1, 4, x[n-3:end], f, x[n])
    w[4,n] = t2 - t1
    
    return
end

