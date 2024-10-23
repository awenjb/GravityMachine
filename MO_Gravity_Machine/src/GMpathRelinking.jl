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


# Fonction de Path Relinking entre deux solutions binaires
function path_relinking(solution1::tSolution{Int64}, solution2::tSolution{Int64}, c1::Array{Int,1}, c2::Array{Int,1}, A,  mode)

	if mode=="N" 
		@info "Mode Naive"
		# Obtient un chemin
		path = naive_path_relinking(solution1, solution2, c1, c2)
		
		var = length(path)
		@info "Génération de $var solutions"
		
		supr = []
		# Vérifier Admissibilité
		for i in eachindex(path)
			if check_admissibility(path[i], A) != true
				push!(supr, i)
			end
		end
		deleteat!(path, supr)

		supr = []
		# Vérifier Dominance(regarder chaque paire de solutions)
		for i in eachindex(path)
			for j in eachindex(path)
				if i != j && check_weak_dominance(path[i], path[j])
					push!(supr, j)
					break
				end
			end
		end
		unique!(supr)
		deleteat!(path, supr)

		var = length(path)-2
		@info "$var nouvelles solutions admissibles non dominées trouvées"

		return path
	end

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
				break
			end
		end
	end
	unique!(supr)
	deleteat!(solutions, supr)
	return solutions
end

#=
c1, c2, A = loadInstance2SPA("didactic5.txt") 

@show size(A)
@show A
nbvar = size(A)[2]

# Génération de deux solutions distinctes
ok0 =  tSolution{Int64}(ones(Int64,nbvar),zeros(Int64,2))
ok1 =  tSolution{Int64}(ones(Int64,nbvar),zeros(Int64,2))

ok0.x = [0,0,1,0,0,1,1,0]
ok1.x = [1, 0, 0, 1, 1, 0, 0, 0]
ok0.y = evaluerSolution(ok0.x, c1, c2)
ok1.y = evaluerSolution(ok1.x, c1, c2)

@show ok0, ok1

path_relinking(ok0, ok1, c1, c2, A, "N")

=#

