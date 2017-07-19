import PyPlot
const plt = PyPlot

function M(p::Int64,u)
    if p == 2
        if u >= 0 && u <= 2
            return 1 - abs(u-1)
        else
            return 0.0
        end
    else
        return M(p-1,u)*u/(p-1) + (p-u)/(p-1)*M(p-1,u-1)
    end
end


function main()
    @time println( M(2,1.1) )
    @time println( M(3,1.1) )
    @time println( M(4,1.1) )

    NptsPlot = 100
    x = Array{Float64}( linspace(0.0,5.0,NptsPlot) )
    y2 = Array{Float64}(NptsPlot)
    y3 = Array{Float64}(NptsPlot)
    y4 = Array{Float64}(NptsPlot)
    y5 = Array{Float64}(NptsPlot)
    for i = 1:NptsPlot
        y2[i] = M(2,x[i])
        y3[i] = M(3,x[i])
        y4[i] = M(4,x[i])
        y5[i] = M(5,x[i])
    end
    plt.clf()
    plt.plot( x, y2, linewidth=2.0, label="M2.png" )
    plt.plot( x, y3, linewidth=2.0, label="M3.png" )
    plt.plot( x, y4, linewidth=2.0, label="M4.png" )
    plt.plot( x, y5, linewidth=2.0, label="M5.png" )
    plt.grid(true)
    plt.legend()
    plt.savefig("M2.png", dpi=300)
end

main()
