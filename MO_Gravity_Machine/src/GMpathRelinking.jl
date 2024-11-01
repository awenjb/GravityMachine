include("GMdatastructures.jl")
# ==============================================================================
# Path relinking between two solution

# Utilise un 0-1 échange pour aller de la solution 1 vers la solution 2
# Retourne un chemin de la solution 1 vers la solution 2
function naive_path_relinking(solution1::tSolution{Int64}, solution2::tSolution{Int64}, c1::Array{Int,1}, c2::Array{Int,1})

	# Solution courrante
	current_solution = deepcopy(solution1)

	# Solutions intermediaires
	intermediate_solution = []
	push!(intermediate_solution, solution1)

	# Indice des éléments différents
	diff_indices = findall(x -> solution1.x[x] != solution2.x[x], 1:length(solution1.x))

	# Explorer le chemin vers la solution guide
	for i in diff_indices

		# Modifier l'indice
		current_solution.x[i] = solution2.x[i]
		
		# Calcul du coût
		if solution2.x[i] == 1 
			current_solution.y[1] += c1[i]
			current_solution.y[2] += c2[i]
		else 
			current_solution.y[1] -= c1[i]
			current_solution.y[2] -= c2[i]
		end

		# Ajout au chemin
		push!(intermediate_solution, deepcopy(current_solution))
	end
	return intermediate_solution
end


# Utilise un 0-1 échange pour aller de la solution 1 vers la solution 2
# Retourne tous les chemins de la solution 1 vers la solution 2 possibles
function brute_path_relinking(solution1::tSolution{Int64}, solution2::tSolution{Int64}, c1::Array{Int,1}, c2::Array{Int,1})

    # Solutions intermédiaires
    intermediate_solutions = []
    push!(intermediate_solutions, deepcopy(solution1))

    # Indice des éléments différents
    diff_indices = findall(x -> solution1.x[x] != solution2.x[x], 1:length(solution1.x))

    # Nombre de solutions intermédiaires
    n_changes = length(diff_indices)

    # Parcourir toutes les combinaisons possibles de changements
    for k in 0:(2^n_changes - 1)
        # Créer une copie de la solution initiale
        current_solution = deepcopy(solution1)

        # Parcourir 
        for j in 1:n_changes
            # Vérifier si la variable j du nombre k est à 1
            if (k >> (j - 1)) & 1 == 1

                index = diff_indices[j]
                current_solution.x[index] = solution2.x[index]

                # Mise à jour des coûts
                if solution2.x[index] == 1
                    current_solution.y[1] += c1[index]
                    current_solution.y[2] += c2[index]
                else
                    current_solution.y[1] -= c1[index]
                    current_solution.y[2] -= c2[index]
                end
            end
        end
        
        # Ajout de la solution intermédiaire au chemin
        push!(intermediate_solutions, deepcopy(current_solution))
    end

	# supprimer doublons
	supr = []
	for i in eachindex(intermediate_solutions)
        for j in (i+1):length(intermediate_solutions)
            if intermediate_solutions[i].x == intermediate_solutions[j].x
                push!(supr, j)
            end
        end
    end
    deleteat!(intermediate_solutions, supr)

    return intermediate_solutions
end	

# Fonction de Path Relinking entre deux solutions binaires
function path_relinking(solution1::tSolution{Int64}, solution2::tSolution{Int64}, c1::Array{Int,1}, c2::Array{Int,1}, A,  mode)

	path::Array{tSolution{Int64}} = []
	
	if mode=="N" 
		@info "Mode Naive"
		# Obtient un chemin
		path = naive_path_relinking(solution1, solution2, c1, c2)
	
	end	
	
	if mode=="B" 
		@info "Mode BruteForce"
		# Obtient un chemin
		path = brute_path_relinking(solution1, solution2, c1, c2)
	
	end	
	var = length(path)
		
	@info "Génération de $var solutions"
		
	supr = []
	# Vérifier Admissibilité
	for i in eachindex(path)
		if isFeasible(path[i], A) != true
			push!(supr, i)
		end
	end
	unique!(supr)
	deleteat!(path, supr)

	path = remove_dominated(path)

	var = length(path)-2
	@info "$var nouvelles solutions admissibles non dominées trouvées"

	return path

end


function evaluerSolution(x::Vector{Int64}, c1::Array{Int,1}, c2::Array{Int,1})
	
    z1 = 0.0; z2 = 0.0
    for i in 1:length(x)
        z1 += x[i] * c1[i]
        z2 += x[i] * c2[i]
    end
    return [round(z1, digits=2), round(z2, digits=2)]
end

# check admissibility for a solution vector of SPA
function check_admissibility(solution::tSolution{Int64}, A::Matrix{Int})

    nbctr, nbvar = size(A)
    coverage = zeros(Int, nbctr)

    for j in 1:nbvar
        if solution.x[j] == 1
            coverage .+= A[:, j]
        end
    end

    for i in 1:nbctr
        if coverage[i] != 1
            return false
        end
    end
    return true
end

# Vérifie si la solution 1 est meilleure que la solution 2
# utilise les définitions du cours de MO
function check_dominance(solution1::tSolution{Int64}, solution2::tSolution{Int64})
	better_or_equal = all(solution1.y .<= solution2.y)
	strictly_better = any(solution1.y .< solution2.y)
	return better_or_equal && strictly_better
end

function check_weak_dominance(solution1::tSolution{Int64}, solution2::tSolution{Int64})
	return all(solution1.y .<= solution2.y)
end

function remove_dominated(solutions::Array{tSolution{Int64}})
	supr = []
	for i in eachindex(solutions)
		for j in eachindex(solutions)
			if i != j && check_weak_dominance(solutions[i], solutions[j])
				push!(supr, j)
			end
		end
	end
	unique!(supr)
	deleteat!(solutions, supr)
	return solutions
end

#=============================
	Biased Path Relinking
=============================#

"""
# Arguments :
- `U` : The upper bounding set of solutions.
- `roundedSols` : The set of non feasible solutions considered.
- `numberOfRoundedSolPerIteration::Int` : The number of rounded solutions to be considered as bias at for each pair of source/target solutions.
- `A` : The constraint matrix of the studied SPA instance.
- `c1::Vector` : The vector of costs of the first objective function.
- `c2::Vector` : The vector of costs of the second objective function.
- `biasPermutation::tSolution` : The permutation of the order the indices are considered during the biasing process.
- `biasAmount::tSolution` : The amount of indices that should be considered in the bias permutation, also maximum amount of bias applied.
- `pathRelinkingPermutation::tSolution` : The permutation of the order the indices are considered during the path relinking process.
"""
function biased_path_relinking(U::Vector{tSolution{Int}}, roundedSols::Vector{tSolution{Int}}, numberOfRoundedSolPerIteration::Int, A, c1 ,c2, biasPermutation::Vector{Int}=collect(1:length(biasSol.x)), biasAmount::Int=length(biasSol.x), pathRelinkingPermutation=collect(1:length(sourceSol.x)))
	U = [update(e,c1,c2) for e in remove_duplicates(U)] # Remove duplicates in U and update the Y values of the Upper bound solutions
	
	# Create pairs
	sourceTargetPairs = [[u1, U[i2]] for (i1,u1) in enumerate(U) for i2 in i1+1:length(U) if u1.x != U[i2].x]

	# Initialize the resulting set of feasible solutions found with the path relinking
	resultFeasibleSet = Vector{tSolution{Int64}}(undef,0)

	for pair in sourceTargetPairs
		sourceSol = pair[1]
		targetSol = pair[2]
		# Sort the rounded sols
		sortedRoundedSols = phase2_biased_path_relinking(roundedSols, sourceSol, targetSol, optimism, c1, c2)
		sortedRoundedSols = sortedRoundedSols[1:numberOfRoundedSolPerIteration]
		
		for biasSol in sortedRoundedSols
			phase3Sols = phase3_biased_path_relinking(sourceSol, targetSol, biasSol, biasPermutation, biasAmount, pathRelinkingPermutation)
			feasibleSols = [s for s in phase3Sols if isFeasible(s, A)]
			vcat([resultFeasibleSet, feasibleSols])
		end
	end
	
	return resultFeasibleSet
end


# Phase 1 : Select the rounded sols to evaluate for the pair
"""
Sorts rounded solutions according to their distance to a theta point for a given source and target solution, and an optimism parameter.

# Arguments :
- `roundedSols` : The set of non feasible solutions considered.
- `sourceSol::tSolution` : The solution from which the PR starts.
- `targetSol::tSolution` : The solution with which the PR ends.
- `optimism::Float64` : A value between 0 and 1 the linear interpolation relies on.
- `c1::Vector` : The vector of costs of the first objective function.
- `c2::Vector` : The vector of costs of the second objective function.
"""
function phase1_biased_path_relinking(roundedSols::Vector{tSolution}, sourceSol::tSolution, targetSol::tSolution, optimism::Float64, c1, c2)
	remove_duplicates(roundedSols) # Ensures there are no duplicates
	θ = theta_point(sourceSol, targetSol, optimism) # Computes the θ point
	sortedRoundedSols = sort(roundedSols, θ, c1, c2) # Sort the roundedSols according to their objective distance to the θ point
	return sortedRoundedSols
end


# Phase 2 : For each rounded sol to evaluate as the bias, apply tri_path_relinking
"""
Creates new solutions found along the path between a source solution and a target solution, the path being biased by another solution.

# Arguments :
- `sourceSol::tSolution` : The solution from which the PR starts.
- `targetSol::tSolution` : The solution with which the PR ends.
- `biasSol::tSolution` : The solution that biases the PR.
- `biasPermutation::tSolution` : The permutation of the order the indices are considered during the biasing process.
- `biasAmount::tSolution` : The amount of indices that should be considered in the bias permutation, also maximum amount of bias applied.
- `pathRelinkingPermutation::tSolution` : The permutation of the order the indices are considered during the path relinking process.
"""
function phase2_biased_path_relinking(sourceSol::tSolution, targetSol::tSolution, biasSol::tSolution, biasPermutation::Vector{Int}=collect(1:length(biasSol.x)), biasAmount::Int=length(biasSol.x), pathRelinkingPermutation=collect(1:length(sourceSol.x)))
	indexes = collect(1:length(sourceSol))
	indexesNotEqualSourceTarget = [i for i in indexes if sourceSol.x[i] != targetSol.x[i]]
	indexesEqualSourceTarget = [i for i in indexes if sourceSol.x[i] == targetSol.x[i]]
	indexesNotEqualBias = [i for i in indexesEqualSourceTarget if sourceSol.x[i] != biasSol.x[i]]

	sort_by_permutation!(indexesNotEqualBias, biasPermutation)
	indexesNotEqualBias = indexesNotEqualBias[1:biasAmount]

	# Creation of new source and target solutions from the selected indexes of the bias solution
	newSource::tSolution = deepcopy(sourceSol)
	newTarget::tSolution = deepcopy(targetSol)
	for i in indexesNotEqualBias
		newSource.x[i] = biasSol.x[i]
		newTarget.x[i] = biasSol.x[i]
	end

	# Actual path relinking method
	sort_by_permutation!(indexesNotEqualSourceTarget, pathRelinkingPermutation)
	newSolutions = Vector{tSolution}(undef,0)
	currentSol = deepcopy(newSource)
	for i in indexesNotEqualSourceTarget
		currentSol = deepcopy(currentSol)
		currentSol.x[i] = newTarget.x[i]
		push!(newSolutions, currentSol)
	end

	return newSolutions
end


#==================================
	Biased Path Relinking UTILS
===================================#

# Function to create inverse permutation lookup
function inverse_permutation(P::Vector{Int})
    n = length(P)
    inv_P = zeros(Int, n)
    for i in 1:n
        inv_P[P[i]] = i
    end
    return inv_P
end

# Main function to sort v according to permutation P
function sort_by_permutation!(v::Vector{Int}, P::Vector{Int})
    # Get inverse permutation for O(1) lookup
    inv_P = inverse_permutation(P)
    # Sort v based on the positions in inv_P
    return sort!(v, by=x -> inv_P[x])
end

"""
From a vector of solutions, remove the duplicates
"""
function remove_duplicates(sols::Vector{tSolution{Int}})
	a = [e.x for e in sols]
	indexes = unique(i -> a[i], eachindex(a))
	return sols[indexes]
end

# Replace check_admissibility by isFeasible in all code
"""
Tests wether a solution 
"""
function isFeasible(sol::tSolution, A::Matrix)
	return all(A*sol.x .== 1)
end

"""
Evaluates the objective space values of a solution.
"""
function evaluate(s::tSolution, c1::Vector, c2::Vector)
	return [s.x'*c1, s.x'*c2]
end

"""
Updates the objective space values of a solution and returns the evaluated solution.
"""
function update(s::tSolution, c1::Vector, c2::Vector)
	s.y = evaluate(s,c1,c2)
	return s
end

"""
Updates the objective space values of a solution and returns the evaluated solution.
"""
function update!(s::tSolution, c1::Vector, c2::Vector)
	s.y = evaluate(s,c1,c2)
end

"""
Returns the ideal point betweeen two evaluated solutions.
"""
function ideal_point(s1::tSolution, s2::tSolution)
	return [minimum([s1.y[i], s2.y[i]]) for i in 1:length(s1.y)]
end

"""
Returns the Nadir point betweeen two evaluated solutions.
"""
function nadir_point(s1::tSolution, s2::tSolution)
	return [maximum([s1.y[i], s2.y[i]]) for i in 1:length(s1.y)]
end

"""
Returns the euclidean distance between the evaluation of a solution in the objective space and a point in the objective space.
"""
function distance_euclidean(s::tSolution, point::Vector)
	return sqrt(sum((s.y[i] - point[i])^2 for i in 1:length(point)))
end

"""
Computes the θ point using linear interpolation and the optimism factor between 0 and 1 (representing how close to the ideal point)
"""
function theta_point(s1::tSolution, s2::tSolution, optimism::Float64)
	ideal = ideal_point(s1,s2)
	nadir = nadir_point(s1,s2)
	return nadir + optimism*(ideal - nadir)
end

"""
From a vector of rounded points solutions, orders them according to a distance to a θ point.
Operates an update on the objective values of the solutions provided in L for safety.
"""
function sort(L::Vector{tSolution{Int}}, thetaPoint::Vector, c1, c2)
	# Makes sure each point is updated
	for s in L
		update!(s,c1,c2)
	end
	# Compute the distance of each point of s to the θ point
	distances = [distance_euclidean(s,thetaPoint) for s in L]
	# Sort the permutation of the distances
	perms = sortperm(distances)
	
	return L[perms]
end

"""
From a vector of rounded points solutions, orders them according to a distance to a θ point.
!!! No check on the objective values of the solutions of L.
"""
function sort(L::Vector{tSolution{Int}}, thetaPoint::Vector)
	# Compute the distance of each point of s to the θ point
	distances = [distance_euclidean(s,thetaPoint) for s in L]
	# Sort the permutation of the distances
	perms = sortperm(distances)
	
	return L[perms]
end




#=

# from src/ folder
c1, c2, A = loadInstance2SPA("didactic5.txt") 

@show size(A)
@show A
nbvar = size(A)[2]

# Génération de deux solutions distinctes
ok0 =  tSolution{Int64}(ones(Int64,nbvar),zeros(Int64,2))
ok1 =  tSolution{Int64}(ones(Int64,nbvar),zeros(Int64,2))

ok0.x = [0, 0, 1, 0, 0, 1, 1, 0]
ok1.x = [1, 0, 0, 0, 1, 0, 1, 0]
ok0.y = evaluerSolution(ok0.x, c1, c2)
ok1.y = evaluerSolution(ok1.x, c1, c2)

@show ok0, ok1

path_relinking(ok0, ok1, c1, c2, A, "B")


=#