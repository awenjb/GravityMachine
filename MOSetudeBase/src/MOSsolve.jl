
# ==============================================================================
# compute S_N for a 2SPA (JuMP/MOA model) for methods and options selected

function solve2SPA(  
    m2SPA::Model,          # a 2SPA model
    solverMIP::Symbol,     # MIP solver to use (:GLPK or :Gurobi)
    methodMOA::Symbol,     # MOA method to use (:EpsilonConstraint)
    varType::Symbol;       # nature of variables (:Bin or :Con)
    nbPoints::Int64=0      # number of points to compute (optional parameter with default value)
    )


    # ---- precondition
    if methodMOA==:EpsilonConstraint && varType==:Con && nbPoints==0
        @assert false "warning: no reasonable to probe all YN with ϵ-constraint and 0≤x≤1"        
    end        


    # ---- Get the dimensions of the 2SPA instance to solve
    nbvar = num_variables(m2SPA)
    nbctr = num_constraints(m2SPA, AffExpr, MOI.EqualTo{Float64})
    start = time()


    # ---- Relax the integrality constraints on variables if varType == :Con
    if varType == :Con
        undo_relax = relax_integrality(m2SPA)
    end


    # ---- Setting the solver
    if solverMIP == :GLPK
        set_optimizer(m2SPA, () -> MOA.Optimizer(GLPK.Optimizer))
    elseif solverMIP == :Gurobi
        set_optimizer(m2SPA, () -> MOA.Optimizer(Gurobi.Optimizer))
    else
        @assert false "error: unavailable MIP solver requested"
    end   

    set_silent(m2SPA)

    if methodMOA == :EpsilonConstraint
        set_optimizer_attribute(m2SPA, MOA.Algorithm(), MOA.EpsilonConstraint())
    else
        @assert false "error: unavailable MOA method requested"
    end

    if nbPoints > 0
        set_optimizer_attribute(m2SPA, MOA.SolutionLimit(), nbPoints)
    end


    # ---- Run the solver
    optimize!(m2SPA)


    # ---- Querying the results

    cardSN = result_count(m2SPA)
    verbose ? println("  cardSN = $cardSN") : nothing

    SN = Array{Number}(undef,2,cardSN)
    sumNbFrac = 0.0
    for i in 1:cardSN

        if varType == :Bin
            SN[1,i] = convert(Int64,value(m2SPA[:obj1]; result = i))
            SN[2,i] = convert(Int64,value(m2SPA[:obj2]; result = i))
            verbose ? @printf("  %3d: z=[%6d,%6d] | ", i, SN[1,i], SN[2,i]) : nothing
        else
            SN[1,i] = value(m2SPA[:obj1]; result = i)
            SN[2,i] = value(m2SPA[:obj2]; result = i)
            verbose ? @printf("  %3d: z=[%9.2f,%9.2f] | ", i, SN[1,i], SN[2,i]) : nothing
        end

    end


    # ---- Restore the integrality constraints on variables if varType == :Con    
    if varType == :Con
        undo_relax()
    end


    elapsedTime = time()-start
    println("  Elapsed time: $(round(elapsedTime,digits=3))s \n\n ")

    return SN
end

