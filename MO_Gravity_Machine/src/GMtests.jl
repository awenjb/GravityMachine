println("INIT  :  INCLUDING MAIN ")

include("GMmain.jl")
include("GMparsers.jl")

function test_instance_pathrelinking(path::String, samplingSize, maxTrial, maxTime)
    files = readdir(path, join=false)
    for file in files
        # Loads an instance from the instances
        c1, c2, A = loadInstance_2SPA_bis(path*"/"*file) 

        # Launch GM on it
        vU = GM_noGraph(c1, c2, A, samplingSize, maxTrial, maxTime)

        # PATH RELINKING
        U = deepcopy(vU)

        for i in 1:(length(vU)-1)

            PR_path = path_relinking(vU[i], vU[i+1], c1, c2, A, "N")
            
            # Remove initial and target solutions from 
            popfirst!(PR_path)
            pop!(PR_path)

            # add new solutions to U
            for sol in PR_path
                push!(U, sol)
            end
        end
    end
end


function GM_noGraph(c1, c2, A, tailleSampling::Int64, maxTrial::Int64, maxTime::Int64)

    @assert tailleSampling>=3 "Erreur : Au moins 3 sont requis"

    @printf("0) instance et parametres \n\n")
    verbose ? println("  instance = $fname | tailleSampling = $tailleSampling | maxTrial = $maxTrial | maxTime = $maxTime\n\n") : nothing

    # chargement de l'instance numerique ---------------------------------------
    c1, c2, A = loadInstance2SPA(fname) # instance numerique de SPA
    nbctr = size(A,1)
    nbvar = size(A,2)
    nbobj = 2

    # structure pour les points qui apparaitront dans l'affichage graphique
    d = tListDisplay([],[], [],[], [],[], [],[], [],[], [],[], [],[])

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("1) calcule les etendues de valeurs sur les 2 objectifs\n\n")

    # calcule la valeur optimale relachee de f1 seule et le point (z1,z2) correspondant
    f1RL, xf1RL = computeLinearRelax2SPA(nbvar, nbctr, A, c1, c2, typemax(Int), 1) # opt fct 1
    minf1RL, maxf2RL = evaluerSolution(xf1RL, c1, c2)

    # calcule la valeur optimale relachee de f2 seule et le point (z1,z2) correspondant
    f2RL, xf2RL = computeLinearRelax2SPA(nbvar, nbctr, A, c1, c2, typemax(Int), 2) # opt fct 2
    maxf1RL, minf2RL = evaluerSolution(xf2RL, c1, c2)

    verbose ? @printf("  f1_min=%8.2f ↔ f1_max=%8.2f (Δ=%.2f) \n",minf1RL, maxf1RL, maxf1RL-minf1RL) : nothing
    verbose ? @printf("  f2_min=%8.2f ↔ f2_max=%8.2f (Δ=%.2f) \n\n",minf2RL, maxf2RL, maxf2RL-minf2RL) : nothing


    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("2) calcule les generateurs par e-contrainte alternant minimiser z1 et z2\n\n")

    nbgen, L = calculGenerateurs(A, c1, c2, tailleSampling, minf1RL, maxf2RL, maxf1RL, minf2RL, d)

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # allocation de memoire pour la structure de donnees -----------------------

    vg = allocateDatastructure(nbgen, nbvar, nbobj)

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("3) place L dans structure et verifie l'admissibilite de chaque generateur\n\n")

    for k=1:nbgen

    verbose ? @printf("  %2d  : [ %8.2f , %8.2f ] ", k, L[k].y[1], L[k].y[2]) : nothing

    # copie de l'ensemble bornant inferieur dans la stru de donnees iterative ---
    ajouterX0!(vg, k, L[k])


    # test d'admissibilite et marquage de la solution le cas echeant -------
    if estAdmissible(vg[k].sRel.x)
    ajouterXtilde!(vg, k, convert.(Int, vg[k].sRel.x), convert.(Int, L[k].y))
    vg[k].sFea   = true
    verbose ? @printf("→ Admissible \n") : nothing
    # archive le point obtenu pour les besoins d'affichage    
    if generateurVisualise == -1 
        # archivage pour tous les generateurs
        push!(d.XFeas,vg[k].sInt.y[1])
        push!(d.YFeas,vg[k].sInt.y[2])
    elseif generateurVisualise == k
        # archivage seulement pour le generateur k
        push!(d.XFeas,vg[k].sInt.y[1])
        push!(d.YFeas,vg[k].sInt.y[2])
    end 
    else
    vg[k].sFea   = false
    verbose ? @printf("→ x          \n") : nothing
    end

    end
    verbose ? println("") : nothing

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # Sortie graphique

    figure("Gravity Machine",figsize=(6.5,5))
    #xlim(25000,45000)
    #ylim(20000,40000)
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")
    PyPlot.title("Cone | 1 rounding | 2-$fname")

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # calcule les directions (λ1,λ2) pour chaque generateur a utiliser lors des projections
    λ1,λ2 = calculerDirections2(L,vg)

    # ==========================================================================

    @printf("4) terraformation generateur par generateur \n\n")

    for k in [i for i in 1:nbgen if !isFeasible(vg,i)]
    temps = time()
    trial = 0
    H =(Vector{Int64})[]

    #perturbSolution30!(vg,k,c1,c2,d)

    # rounding solution : met a jour sInt dans vg --------------------------
    #roundingSolution!(vg,k,c1,c2,d)  # un cone
    #roundingSolutionnew24!(vg,k,c1,c2,d) # deux cones
    roundingSolutionNew23!(vg,k,c1,c2,d) # un cone et LS sur generateur

    push!(H,[vg[k].sInt.y[1],vg[k].sInt.y[2]])
    println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4))

    while !(t1=isFeasible(vg,k)) && !(t2=isFinished(trial, maxTrial)) && !(t3=isTimeout(temps, maxTime))

    trial+=1

    # projecting solution : met a jour sPrj, sInt, sFea dans vg --------
    projectingSolution!(vg,k,A,c1,c2,λ1,λ2,d)
    println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4))

    if !isFeasible(vg,k)

        # rounding solution : met a jour sInt dans vg --------------------------
        #roundingSolution!(vg,k,c1,c2,d)
        #roundingSolutionnew24!(vg,k,c1,c2,d)
        roundingSolutionNew23!(vg,k,c1,c2,d)
        println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4))

        # test detection cycle sur solutions entieres ------------------
        cycle = [vg[k].sInt.y[1],vg[k].sInt.y[2]] in H
        if (cycle == true)
            println("CYCLE!!!!!!!!!!!!!!!")
            # perturb solution
            perturbSolution30!(vg,k,c1,c2,d)
        end
        push!(H,[vg[k].sInt.y[1],vg[k].sInt.y[2]])

    end
    end
    if t1
    println("   feasible \n")
    elseif t2
    println("   maxTrial \n")
    elseif t3
    println("   maxTime \n")
    end


    end

    println("");


    # ==========================================================================

    @printf("5) Extraction des resultats\n\n")


    for k=1:nbgen
    verbose ? @printf("  %2d  : [ %8.2f , %8.2f ] ", k, vg[k].sInt.y[1],vg[k].sInt.y[2]) : nothing
    # test d'admissibilite et marquage de la solution le cas echeant -------
    if vg[k].sFea
    verbose ? @printf("→ Admissible \n") : nothing
    else
    verbose ? @printf("→ x          \n") : nothing
    end
    end

    # allocation de memoire pour les ensembles bornants ------------------------
    U = Vector{tSolution{Int64}}(undef,nbgen)
    for j = 1:nbgen
    U[j] = tSolution{Int64}(zeros(Int64,nbvar),zeros(Int64,nbobj))
    end
    #--> TODO : stocker l'EBP dans U proprement


    # ==========================================================================
    # Stockage à part des valeurs des variables  des points réalisables de vg pour path relinking

    # vU ensemble bornant supérieur (solutions non dominées)


    vU = Vector{tSolution{Int64}}(undef,nbgen)
    for j = 1:nbgen
    vU[j] = tSolution{Int64}(zeros(Int64,nbvar),zeros(Int64,nbobj))
    end

    vToSupr = zeros(Int64, 0)

    for k=1:nbgen
    # test d'admissibilite et marquage de la solution le cas echeant -------
    if vg[k].sFea
    verbose ? @printf("→ Admissible \n") : nothing
    vU[k].x = vg[k].sInt.x
    vU[k].y = vg[k].sInt.y
    else
    verbose ? @printf("→ x          \n") : nothing
    push!(vToSupr, k)
    end
    end


    # Detection des doublons
    unique!(vU)
    for i=1:size(vU)[1] 
    for j=(i+1):size(vU)[1] 
    if vU[i].x == vU[j].x
        push!(vToSupr, j)
    end
    end
    end

    # Suppression des éléments non réalisables / doublons"
    deleteat!(vU, vToSupr)

    vU = remove_dominated(vU)

    return vU

    #=
    # Pretty print U 
    print("-------- U \n")
    for sol in U
    print(sol.x, "\n")
    end
    =#


end


c1, c2, A = loadInstance2SPA("didactic5.txt") 
# Launch GM on it
# Launch pathrelinking 