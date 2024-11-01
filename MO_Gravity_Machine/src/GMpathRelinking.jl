# ==============================================================================
# Path relinking entre deux solutions
# Retourne un chemin entre les deux solutions (comprises)


# ==============================================================================
# Path relinking simple
# Utilise un 0-1 échange pour aller de la solution 1 vers la solution 2 créant un chemin aléatoire
# Ne garde que les solutions non dominées admissibles
function simple_path_relinking(solution1::tSolution{Int64}, solution2::tSolution{Int64}, c1::Array{Int,1}, c2::Array{Int,1}, A)

    best_solution = solution1
    current_solution = deepcopy(solution1)
    
	# solutions intermediaires (chemin à retourner)
	intermediate_solution = []

	# ajout de la solution initiale pour commencer le chemin
	push!(intermediate_solution, deepcopy(solution1))

    while current_solution.x != solution2.x
        # trouver les indices qui diffèrent entre la solution1 (initiale) et solution2 (cible)
        diff_indices = findall(x -> solution1.x[x] != solution2.x[x], 1:length(solution1.x))
        
        # si différences, appliquer un changement pour se rapprocher de solution2
        if !isempty(diff_indices)

            # choisir un indice au hasard et le remplacer par sa valeur dans solution2
            idx_change = rand(diff_indices)
			current_solution.x[idx_change] = solution2.x[idx_change]
			
			# calcul du coût (à opti)
			current_solution.y = evaluerSolution(current_solution.x, c1, c2)

            # mettre à jour la meilleure solution si amélioration
            if check_weak_dominance(current_solution, best_solution)
				if check_admissibility(current_solution, A)
					#println("est admissible et non dominé !")
					best_solution = deepcopy(current_solution)
					#@show best_solution

					# si la solution est intéressante, update solutions intermédiaires
					push!(intermediate_solution, deepcopy(current_solution))
				end
            end
			# enlever l'indice de la liste des indices différents
			filter!(x -> x != idx_change, diff_indices)

        else
            break  # plus rien à modifier, on s'arrête (boucle for de taille diff_indices ?)
        end
    end

	# ajout de la solution cible pour compléter le chemin
	push!(intermediate_solution, deepcopy(solution2))

	#@show intermediate_solution
    return intermediate_solution
end


# ==============================================================================
# Path relinking naif
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
		current_solution.y = evaluerSolution(current_solution.x, c1, c2)

		# Ajout au chemin
		push!(intermediate_solution, deepcopy(current_solution))
	end
	return intermediate_solution
end

# ==============================================================================
# Path relinking énumératif
# Utilise un 0-1 échange
# Retourne toutes les solutions des chemins de la solution 1 vers la solution 2
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
                current_solution.y = evaluerSolution(current_solution.x, c1, c2)
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
	unique!(supr)
    deleteat!(intermediate_solutions, supr)

    return intermediate_solutions
end	

# ==============================================================================
# Path relinking avec résolution de deux sous problèmes
# faire un sous pb et le résoudre pour compléter les solutons intermédiaires.
# -> semble retourner soit la solution cible/ soit la solution initiale
function heuristic1_path_relinking(solution1::tSolution{Int64}, solution2::tSolution{Int64}, c1::Array{Int,1}, c2::Array{Int,1}, A)

    # Solutions intermédiaires
    intermediate_solutions = []
    push!(intermediate_solutions, deepcopy(solution1))

    # Indice des éléments différents
    diff_indices = findall(x -> solution1.x[x] != solution2.x[x], 1:length(solution1.x))

    # Solution courrante
    current_solution1 = deepcopy(solution1)
	current_solution2 = deepcopy(solution1)

	# resoudre sous problem de SPA avec les indices de diff_indices

	supr::Vector{Int64} = []
	for (i, line) in enumerate(eachrow(A))
		if !line_with_1(line, diff_indices)
			push!(supr, i)
		end
	end

	# réduire la matrice avec seulement les lignes ou des variables de diff_indices apparaissent
    A_reduce = A[setdiff(1:size(A, 1), supr), :]
	# réduire la matrice avec seulement les colonnes ou des variables de diff_indices apparaissent
	A_reduce = A_reduce[:, setdiff(1:size(A, 2), setdiff(1:size(A, 2), diff_indices))]
	c1_reduce = deleteat!(copy(c1), setdiff(1:size(A, 2), diff_indices))
	c2_reduce = deleteat!(copy(c2), setdiff(1:size(A, 2), diff_indices))

	# résoudre le pb soit selon z1, soit selon z2
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
	
	push!(intermediate_solutions, deepcopy(solution2))

	# supprimer doublons
	supr = []
	for i in eachindex(intermediate_solutions)
		for j in (i+1):length(intermediate_solutions)
			if intermediate_solutions[i].x == intermediate_solutions[j].x
				push!(supr, j)
			end
		end
	end
	unique!(supr)
	deleteat!(intermediate_solutions, supr)

    return intermediate_solutions

end	


# ==============================================================================
# Path relinking avec résolution de deux sous problèmes
# Pour la 1ere variable en conflit entre s1 et s2, fixer s1 à la valeure de s2,
# résoudre le sous problème résultant avec relax linéaire puis arrondir pour obtenir une nouvelle solution
# faire ça pour chaque variable ???
function heuristic2_path_relinking(solution1::tSolution{Int64}, solution2::tSolution{Int64}, c1::Array{Int,1}, c2::Array{Int,1}, A)

	# Solution courrante
	current_solution = deepcopy(solution1)

	# Solutions intermediaires
	intermediate_solutions = []
	push!(intermediate_solutions, solution1)

	# Indice des éléments différents
	diff_indices = findall(x -> solution1.x[x] != solution2.x[x], 1:length(solution1.x))

	# Explorer le chemin vers la solution guide
	# pour chaque indice, résoudre un pb relaché

	# récupérer les variables à 1 et 0 communes pour les fixer
	x_s1 = copy(solution1.x)

	v1 = findall(x -> x == 1, x_s1)
	v0 = findall(x -> x == 0, x_s1)
	filter!(x -> !(x in diff_indices), v0)
	filter!(x -> !(x in diff_indices), v1)

	nbvar = size(A,2)
	nbctr = size(A,1)

	n_changes = length(diff_indices)

	for i in 1:n_changes

		println("-----------------------------------")
		obj = 2

		# poser le sous problème
		model = Model(GLPK.Optimizer)
		@variable(model, 0.0 <= x[1:nbvar] <= 1.0 ) # variable reel
		@constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
		# fixer les variables 
		@constraint(model, [i=1:length(v1)], x[v1[i]] == 1)
		@constraint(model, [i=1:length(v0)], x[v0[i]] == 0)

		if obj == 1
		  @objective(model, Min, sum((c1[i])*x[i] for i in 1:nbvar))
		else
		  @objective(model, Min, sum((c2[i])*x[i] for i in 1:nbvar))
		end
		optimize!(model)

		# Récupérer résultats
		current_solution.x = value.(x)  #TODO gérer les approximations et fractions dans .x

		@info "Solutions du soveurs"

		# Récupérer la première valeur la plus proche de 0 ou 1
		check::Vector{Float64} = current_solution.x[diff_indices]
		@show check

		new_index = closer_to_0_1(check)
		# donc diff_indices[new_index] donne l'indice de l'élément à fixer dans x
		@show diff_indices
		@show new_index
		@show diff_indices[new_index]

		@info "Nouvelle solution intermédiaire"

		#var = current_solution.x[diff_indices[new_index]]
		# Rendre la solution entière pour obtenir une solution intermédiaire
		for i in diff_indices
			current_solution.x[i] = round(current_solution.x[i])

			# ajout d'aléatoire ?
			
			#=
			if i != diff_indices[new_index]
				println("aléatoire")
				t = rand([0, 1])
				current_solution.x[i] = t
			end
			=#
		end
		
		current_solution.y = evaluerSolution(current_solution.x, c1, c2)

		@show current_solution.x
		@show current_solution.y

		# Ajouter new_index dans les variables fixées à 0 ou 1
		if current_solution.x[diff_indices[new_index]] == 1
			push!(v1, diff_indices[new_index])
		else
			push!(v0, diff_indices[new_index])
		end
		
		@info "Nouveaux diff indices"

		@show v0,v1
		# supprimer new_index de diff_indices
		#filter!(x -> x != new_index, diff_indices)
		deleteat!(diff_indices, new_index)

		@show diff_indices

		# Ajout au chemin
		push!(intermediate_solutions, deepcopy(current_solution))

	end
	@info "sorti"
	
	push!(intermediate_solutions, solution2)
	@show intermediate_solutions
	# supprimer doublons
	supr = []
	for i in eachindex(intermediate_solutions)
		for j in (i+1):length(intermediate_solutions)
			if intermediate_solutions[i].x == intermediate_solutions[j].x
				push!(supr, j)
			end
		end
	end
	@info "sorti2"
	@show supr
	unique!(supr)
	deleteat!(intermediate_solutions, supr)
	@info "sorti3"
	
	@show intermediate_solutions
	return intermediate_solutions

end	


# ==============================================================================
# Path relinking avec résolution d'un sous problème
# fixe les variables communes des deux solutions
# ajouter des contraintes pour éviter de retomber sur les solutions initiales et cibles
# sens d'optimisation = somme pondérée en fonction des solutions initiales et cibles
# résoudre le sous problème résultant avec relax linéaire puis arrondir pour obtenir une nouvelle solution
# TODO ajouter un ranking pour obtenir plus de solutions intermédiares
function heuristic4_path_relinking(solution1::tSolution{Int64}, solution2::tSolution{Int64}, c1::Array{Int,1}, c2::Array{Int,1}, A)

	# Solution courrante
	current_solution = deepcopy(solution1)

	# Solutions intermediaires
	intermediate_solutions = []
	push!(intermediate_solutions, solution1)

	# Indice des éléments différents
	diff_indices = findall(x -> solution1.x[x] != solution2.x[x], 1:length(solution1.x))

	# récupérer les variables à 1 et 0 communes pour les fixer
	x_s1 = copy(solution1.x)

	v1 = findall(x -> x == 1, x_s1)
	v0 = findall(x -> x == 0, x_s1)
	filter!(x -> !(x in diff_indices), v0)
	filter!(x -> !(x in diff_indices), v1)

	nbvar = size(A,2)
	nbctr = size(A,1)

	λ = (solution1.y[2] - solution2.y[2], solution2.y[1] - solution1.y[1])  

	# poser le sous problème
	model = Model(GLPK.Optimizer)
	@variable(model, 0.0 <= x[1:nbvar] <= 1.0 ) # variable reel
	@constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)

	# fixer les variables 
	@constraint(model, [i=1:length(v1)], x[v1[i]] == 1)
	@constraint(model, [i=1:length(v0)], x[v0[i]] == 0)

	# fixer une epsilon sur chaque objectif
	@constraint(model, sum((c1[i])*x[i] for i in 1:nbvar) <= solution2.y[1]-1) #objectif 1
	@constraint(model, sum((c2[i])*x[i] for i in 1:nbvar) <= solution1.y[2]-1) #objectif 2

	@objective(model, Min, λ[1]*sum((c1[i])*x[i] for i in 1:nbvar) +  λ[2]*sum((c2[i])*x[i] for i in 1:nbvar) )

	optimize!(model)

	# Récupérer résultats
	relax_x = value.(x)
	relax_x_int::Vector{Int64} = []
	for i in eachindex(relax_x)
		push!(relax_x_int, Int64(round(relax_x[i])))
	end

	current_solution = deepcopy(solution1)
	current_solution.x = relax_x_int
	current_solution.y = evaluerSolution(current_solution.x, c1, c2)

	push!(intermediate_solutions, current_solution)
	push!(intermediate_solutions, solution2)

	# supprimer doublons
	supr = []
	for i in eachindex(intermediate_solutions)
		for j in (i+1):length(intermediate_solutions)
			if intermediate_solutions[i].x == intermediate_solutions[j].x
				push!(supr, j)
			end
		end
	end
	unique!(supr)
	deleteat!(intermediate_solutions, supr)

	return intermediate_solutions

end	


# Fonction de Path Relinking entre deux solutions binaires
function path_relinking(solution1::tSolution{Int64}, solution2::tSolution{Int64}, c1::Array{Int,1}, c2::Array{Int,1}, A,  mode)

	path::Array{tSolution{Int64}} = []
	
	if mode=="S" 
		println("Mode Simple")
		path = simple_path_relinking(solution1, solution2, c1, c2, A)
	end	

	if mode=="N" 
		println("Mode Naive")
		path = naive_path_relinking(solution1, solution2, c1, c2)
	end	
	
	if mode=="B" 
		@info "Mode BruteForce"
		path = brute_path_relinking(solution1, solution2, c1, c2)
	end	

	if mode=="H1" 
		@info "Mode Heuristique 1"
		path = heuristic1_path_relinking(solution1, solution2, c1, c2, A)
	end	

	if mode=="H2" 
		@info "Mode Heuristique 2"
		path = heuristic2_path_relinking(solution1, solution2, c1, c2, A)
	end	

	if mode=="H4" 
		@info "Mode Heuristique 4"
		path = heuristic4_path_relinking(solution1, solution2, c1, c2, A)
	end	

	var = length(path)
		
	println("Chemin de taille $var")
		
	supr = []
	# Vérifier Admissibilité
	for i in eachindex(path)
		if isFeasible(path[i], A) != true
			push!(supr, i)
		end
	end
	unique!(supr)
	deleteat!(path, supr)

	# Vérifier Dominance
	path = remove_dominated(path)

	#@show path

	# Solutions initiales et cibles comprises dans le chemin -> -2
	var = length(path)-2
	println("$var nouvelles solutions admissibles non dominées trouvées")

	return path
end

# ==============================================================================
# FONCTIONS UTILES
# ==============================================================================

function evaluerSolution(x::Vector{Int64}, c1::Array{Int,1}, c2::Array{Int,1})
	
    z1 = 0.0; z2 = 0.0
    for i in 1:length(x)
        z1 += x[i] * c1[i]
        z2 += x[i] * c2[i]
    end
    return [round(z1, digits=2), round(z2, digits=2)]
end

# Vérifie l'admissibilité d'une solution pour le problème
# Produit matrice vecteur plus efficace ?
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

function closer_to_0_1(vect::Vector{Float64})
    # Calculer les distances par rapport à 0 et 1
    argmin([min(abs(v - 0), abs(v - 1)) for v in vect])
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
	sort!(supr)
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


# 14 48 
#ok0.x = [0, 0, 1, 0, 0, 1, 1, 0]
# 18 16
ok0.x = [1, 0, 0, 0, 1, 0, 1, 0]
# 20 7
#ok1.x = [1, 0, 0, 1, 1, 0, 0, 0]
#21 4
ok1.x = [0, 0, 0, 0, 1, 0, 0, 1]


#ok0.x = [0, 0, 1, 0, 0, 1, 1, 0]
#ok1.x = [1, 0, 0, 1, 1, 0, 0, 0]


ok0.y = evaluerSolution(ok0.x, c1, c2)
ok1.y = evaluerSolution(ok1.x, c1, c2)

@show ok0, ok1

path_relinking(ok0, ok1, c1, c2, A, "S")
=#