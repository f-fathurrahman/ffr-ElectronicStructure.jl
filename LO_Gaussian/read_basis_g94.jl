module g94basis

struct BasisSetParams
    Nsyms::Int64
    syms::Array{String,1} # syms[1:Nsyms]
    Ncontr::Array{Int64,1} # Ncontr[1:Nsyms]
    expns::Array{Array{Float64,1},1} # expns[1:Nsyms][1:Ncontr]
    coeffs::Array{Array{Float64,1},1} # coeffs[1:Nsyms][1:Ncontr]
end

function my_read_lines(filename)
    f = open(filename, "r")
    lines = ""
    while !eof(f)
        lines = lines*readline(f)*"\n"
    end
    close(f)
    slines = split(lines, "****")
    return slines
end

function read_basis_set(filename)
    slines = my_read_lines(filename)
    #
    BasisSetData = Dict{String,BasisSetParams}()
    #
    for l in slines
        if ! ( (l == "\n") || (l == "") )
            basis_params = process_lines(l)
            push!( BasisSetData, basis_params )
        else
            continue    
        end
    end
    return BasisSetData
end

function calc_Nsyms(lines)
    ls = split(lines, "\n")
    # start from line 2
    il = 2
    Nsyms = 0
    while il < length(ls)-1
        il = il + 1
        #
        sym_read = split(ls[il])[1]
        #
        ncontr = parse( Int64, split(ls[il])[2] )
        #
        if sym_read == "SP"
            Nsyms = Nsyms + 2
            il = il + ncontr
        elseif (sym_read == "S") || (sym_read == "P") ||
               (sym_read == "D") || (sym_read == "F")
            Nsyms = Nsyms + 1
            il = il + ncontr
        else
            println("Unknown sym_read: ", sym_read)
            exit()
        end
    end
    return Nsyms
end

function process_lines(lines)

    Nsyms = calc_Nsyms(lines)

    ls = split(lines, "\n")
    
    # start from line 2
    atsymb = split( ls[2] )[1]

    il = 2
    
    Ncontr = Array{Int64}(undef,Nsyms)
    syms = Array{String}(undef,Nsyms)

    expns = Array{Array{Float64,1},1}(undef,Nsyms)
    coeffs = Array{Array{Float64,1},1}(undef,Nsyms)    
    
    isym = 0

    while il < length(ls)-1
        
        il = il + 1
                        
        sym_read = split(ls[il])[1]
        
        ncontr = parse( Int64, split(ls[il])[2] )

        if sym_read == "SP"

            expn = zeros(ncontr)
            coeff_s = zeros(ncontr)
            coeff_p = zeros(ncontr)
            
            for i = 1:ncontr
                #
                il = il + 1 # don't forget to increament index for line
                #
                expn[i] = parse( Float64, split(ls[il])[1] )
                coeff_s[i] = parse( Float64, split(ls[il])[2] )
                coeff_p[i] = parse( Float64, split(ls[il])[3] )
            end
            
            isym = isym + 1 # first we add S
            syms[isym] = "S"
            Ncontr[isym] = ncontr
            expns[isym] = expn[:]
            coeffs[isym] = coeff_s[:]

            isym = isym + 1 # then we add P
            syms[isym] = "P"
            Ncontr[isym] = ncontr
            expns[isym] = expn[:]
            coeffs[isym] = coeff_p[:]  # use the coeff for p

        else 
            # other than sym_read == "SP"

            expn = zeros(ncontr)
            coeff = zeros(ncontr)
            
            for i = 1:ncontr
                #
                il = il + 1 # don't forget to increament index for line
                #
                expn[i] = parse( Float64, split(ls[il])[1] )
                coeff[i] = parse( Float64, split(ls[il])[2] )
            end
            
            isym = isym + 1
            syms[isym] = sym_read
            Ncontr[isym] = ncontr
            expns[isym] = expn[:]
            coeffs[isym] = coeff[:]
        end
    end
    return (atsymb => BasisSetParams(Nsyms, syms, Ncontr, expns, coeffs))
end

end # module g94basis

#import .g94basis
#println( g94basis.read_basis_set("g94_sto-3g") )