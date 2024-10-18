include("GMdatastructures.jl")
# ==============================================================================
# Path relinking between two solution

#
function generate_intermediate_solution(solution1::Vector{Int}, solution2::Vector{Int}, c1::Array{Int,1}, c2::Array{Int,1}, k::Int)

    intermediate_solution = copy(solution1)
    
	# Liste des indices où les deux solutions sont différentes
	diff_indices = findall(x -> solution1[x] != solution2[x], 1:length(solution1))
    
	# Modifie au moins k variables de solution1 pour les rapprocher de solution2
	for i in 1:min(k, length(diff_indices))
    	idx = diff_indices[i]
    	intermediate_solution[idx] = solution2[idx]
	end
    
	return intermediate_solution
end


# Fonction de Path Relinking entre deux solutions binaires
function path_relinking(solution1, solution2, c1::Array{Int,1}, c2::Array{Int,1})

    best_solution = solution1.x
	best_score = solution1.y
    
	# Distance entre les deux solutions
	num_diff = sum(solution1.x .!= solution2.x)
    
	# Génération des solutions intermédiaires en modifiant 1, 2, ..., num_diff variables
	for k in 1:num_diff
    	intermediate_solution = generate_intermediate_solution(solution1.x, solution2.x, k)

        # Vérifier admissibilité
        # Vérifier dominance
	end
    
	return best_solution, best_score
end


# test un deux test
ok0 =  tSolution{Int64}(zeros(Int64,5),zeros(Int64,2))
ok1 =  tSolution{Int64}(ones(Int64,5),ones(Int64,2))

path_relinking(ok0, ok1)