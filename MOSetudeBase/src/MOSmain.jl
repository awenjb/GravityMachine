# ==============================================================================
# ==============================================================================

println("-"^80)
println("""\nMOS2024 \n""")

global verbose    = true
global graphic    = true
global exact      = false
global experiment = false

print("  verbose.....: "); verbose    ? println("yes") : println("no") 
print("  graphics....: "); graphic    ? println("yes") : println("no") 
print("  exact.......: "); exact      ? println("yes") : println("no") 
print("  experiment..: "); experiment ? println("yes") : println("no") 

println("\n\n-) Activation des packages necessaires\n")
using JuMP, GLPK
import MultiObjectiveAlgorithms as MOA
using Printf, PyPlot
println("  Fait\n")

include("MOSfiles.jl")
include("MOSjump.jl")
include("MOSsolve.jl")


# ==============================================================================
# point d'entree principal

function main(fname::String)


    # --------------------------------------------------------------------------
    # Charge une instance numerique   

    println("\n0) instance et dimensions \n")
    println("  instance = $fname") 
   
    m2SPA = load2SPA(fname)
    println("  nbvar    = ",num_variables(m2SPA))
    println("  nbctr    = ",num_constraints(m2SPA, AffExpr, MOI.EqualTo{Float64}),"\n\n")  



    # --------------------------------------------------------------------------
    # Calcul de L pour Y_N (cad les generateurs) 

    println("1) calcule LB(Y_N) avec methode ϵ-constraint et 0≤x≤1")
    nbProbe = 16
    LBE = solve2SPA(m2SPA, :GLPK, :EpsilonConstraint, :Con, nbPoints=nbProbe)
 
    

    # --------------------------------------------------------------------------
    # Sortie graphique

    if graphic

        fig1 = figure("Objective Space Y",figsize=(6.5,5))
        title("single ϵ-constraint | max " * string(nbProbe) * " samples")
        xlabel(L"z^1(x)")
        ylabel(L"z^2(x)")

        # LBE : generateurs obtenus avec une ϵ-constraint simple avec un pas predefini (-> echantillonnage)
        scatter(LBE[1,:], LBE[2,:], c="red", marker="x", s=80, label=L"$LB$ eps") 
        # a faire : ϵ-constraint simple sur l'autre objectif
        legend() 

    end

    return nothing

end


# ==============================================================================
# Run all algorithms over all instances

function numericalExperiment(target)

    global graphic = false
    global exact   = false # some instances are very long to be solved with glpk
    
    fnames = getfname(target)
    for instance = 1:length(fnames)
        start = time()
        main( string(target,"/",fnames[instance]) )    
        elapsedTime = time()-start
        println("Elapsed time for $(fnames[instance]) : $elapsedTime (s)")
    end

    return nothing
end


# ==============================================================================

target = "../SPA/instances"

if experiment
    numericalExperiment(target)
else
    #@time main(target*"/bio"*"sppaa02.txt")
    #@time main(target*"/bio"*"sppnw03.txt")
    #@time main(target*"/bio"*"sppnw04.txt")
    #@time main(target*"/bio"*"sppnw10.txt")
    #@time main(target*"/bio"*"sppnw20.txt")
    #@time main(target*"/bio"*"sppnw25.txt")
    #@time main(target*"/bio"*"didactic3.txt")
    @time main(target*"/bio"*"didactic5.txt")
    #@time main(target*"/bio"*"sppnw29.txt")
    #@time main(target*"/bio"*"sppnw19.txt")
    #@time main(target*"/bio"*"sppnw40.txt")
end

nothing
