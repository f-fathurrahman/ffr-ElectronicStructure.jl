function integrate_trapz_7(f, dr)
    N = size(dr, 1)
    g = f .* dr
    s = (  36799 * ( g[1] + g[N  ] ) +
          176648 * ( g[2] + g[N-1] ) +
           54851 * ( g[3] + g[N-2] ) +
          177984 * ( g[4] + g[N-3] ) +
           89437 * ( g[5] + g[N-4] ) +
          130936 * ( g[6] + g[N-5] ) +
          119585 * ( g[7] + g[N-6] )
    ) / 120960
    for i in 8:N-7
        s = s + g[i]
    end
    return s
end
