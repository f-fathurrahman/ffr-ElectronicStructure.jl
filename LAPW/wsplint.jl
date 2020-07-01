function wsplint!(n, x, w)

    if n <= 9
        for i in 1:n
          f[:] .= 0.0
          f[i] = 1.0
          w[i] = splint(n,x,f)
        end
        return
    end

    f = zeros(9)
    for i in 1:4
        f[:] .= 0.0
        f[i] = 1.0
        w[i] = splint(9, x, f)
    end
    
    f[:] .= 0.0
    f[5] = 1.0
    for i in 5:n-4
      @views w[i] = splint(9, x[i-4:end], f)
    end
    
    for i in 1:4
      f[:] .= 0.0
      f[i+5] = 1.0
      @views w[n-4+i] = splint(9, x[n-8:end], f)
    end
    return
end

