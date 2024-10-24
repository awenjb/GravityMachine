include("GMdatastructures.jl")
include("GMjumpModels.jl")
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

# ça fonctionne pas
function heuristic_path_relinking(solution1::tSolution{Int64}, solution2::tSolution{Int64}, c1::Array{Int,1}, c2::Array{Int,1})

    # Solutions intermédiaires
    intermediate_solutions = []
    push!(intermediate_solutions, deepcopy(solution1))

    # Indice des éléments différents
    diff_indices = findall(x -> solution1.x[x] != solution2.x[x], 1:length(solution1.x))

    # Solution courrante
    current_solution1 = deepcopy(solution1)
	current_solution2 = deepcopy(solution1)

	# TODO 
	# resoudre sous problem de SPA avec les indices de diff_indices

	supr::Vector{Int64} = []
	for (i, line) in enumerate(eachrow(A))
		if !ligne_a_1_dans_indices(line, diff_indices)
			push!(supr, i)
		end
	end

	#indices_filtres = setdiff(1:size(A, 2), diff_indices)

	# réduire la matrice avec seulement les lignes ou des variables de diff_indices apparaissent
    A_reduce = A[setdiff(1:size(A, 1), supr), :]
	# réduire la matrice avec seulement les colonnes ou des variables de diff_indices apparaissent
	A_reduce = A_reduce[:, setdiff(1:size(A, 2), setdiff(1:size(A, 2), diff_indices))]

	c1_reduce = deleteat!(copy(c1), setdiff(1:size(A, 2), diff_indices))
	c2_reduce = deleteat!(copy(c2), setdiff(1:size(A, 2), diff_indices))

	# résoudre le pb 
	zRL1, xRL1 = computeLinearRelax2SPA(size(A_reduce, 2), size(A_reduce, 1), A_reduce, c1_reduce, c2_reduce, typemax(Int), 1) 
	zRL2, xRL2 = computeLinearRelax2SPA(size(A_reduce, 2), size(A_reduce, 1), A_reduce, c1_reduce, c2_reduce, typemax(Int), 2) 
	
	# Arrondir solutions à l'entier le plus proche

	for i in 1:size(A_reduce, 2)
		xRL1[i] = round(xRL1[i])
		xRL2[i] = round(xRL2[i])
	end

	# Créer solution intermédiaire avec les morceaux de solutions trouvés
	k=1
	for i in diff_indices
		current_solution1.x[i] = xRL1[k]
		current_solution2.x[i] = xRL2[k]
		k += 1
	end

	# Evaluer solutions
	current_solution1.y = evaluerSolution(current_solution1.x, c1, c2)
	current_solution2.y = evaluerSolution(current_solution2.x, c1, c2)
	
	push!(intermediate_solutions, current_solution1)
	push!(intermediate_solutions, current_solution2)
	
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

	if mode=="H1" 
		@info "Mode Heuristique 1"
		# Obtient un chemin
		path = heuristic_path_relinking(solution1, solution2, c1, c2)
	end	


	var = length(path)
		
	@info "Génération de $var solutions"
		
	supr = []
	# Vérifier Admissibilité
	for i in eachindex(path)
		if check_admissibility(path[i], A) != true
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

function line_with_1(line, index)
    any(line[idx] == 1 for idx in index)
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

#=
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

path_relinking(ok0, ok1, c1, c2, A, "H1")
=#