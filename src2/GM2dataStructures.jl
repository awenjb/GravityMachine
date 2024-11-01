"""
    Encodes a solution from the feasible space of solutions X and its result
    on the feasible objective space Y.
"""
mutable struct tSolution{T}
    x::Vector{T} # Value in the solution space
    y::Vector{T} # Value in the objective space
end

"""
    Encodes the information of a generator solution.
"""
mutable struct tGenerator
    sRelaxed::tSolution{Float64} # 

    sInteger::tSolution{Int64}
    sIntFeasible::Bool # States if the integer solution is feasible

    sProjected::tSolution{Float64}
end