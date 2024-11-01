
# ==============================================================================
# Create (set) a 2SPA model (bi-objective partionning problem) with JuMP/MOA 
# from a instance file,  with all variables :Bin or :Con

function set2SPA(C::Array{Int,2}, A::Array{Int,2}, varType::Symbol)

    nbctr, nbvar = size(A)
    m2SPA = Model()
    if varType == :Bin
        @variable(m2SPA, x[1:nbvar], Bin)
    elseif varType == :Con
        @variable(m2SPA, 0.0 <= x[1:nbvar] <= 1.0 )
    else
        @assert false "error: unknow type to set for the variables"        
    end
    @expression(m2SPA, obj1, sum((C[1,i])*x[i] for i in 1:nbvar))
    @expression(m2SPA, obj2, sum((C[2,i])*x[i] for i in 1:nbvar))
    @objective(m2SPA, Min, [obj1, obj2])
    @constraint(m2SPA, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1) 

    return m2SPA
end



# ==============================================================================
# Create (load) a 2SPA model (bi-objective partionning problem) with JuMP/MOA 
# from a instance file,  with all variables set to :Bin

function load2SPA(fname::String)

    C, A = parse2SPA(fname)
    m2SPA = set2SPA(C, A, :Bin)

    return m2SPA
end
