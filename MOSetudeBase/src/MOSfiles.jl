# ==============================================================================
# ==============================================================================
# collect the un-hidden filenames available in a given directory

function getfname(target)
    # target : string := chemin + nom du repertoire ou se trouve les instances

    # positionne le currentdirectory dans le repertoire cible
    #cd(joinpath(homedir(),target))
    savePwd = pwd()
    newPwd = pwd() * "/" * target
    cd(newPwd)
    # retourne le repertoire courant
    println("pwd = ", pwd())

    # recupere tous les fichiers se trouvant dans le repertoire data
    allfiles = readdir()

    # vecteur booleen qui marque les noms de fichiers valides
    flag = trues(size(allfiles))

    k=1
    for f in allfiles
        # traite chaque fichier du repertoire
        if f[1] != '.'
            # pas un fichier cache => conserver
            println("fname = ", f)
        else
            # fichier cache => supprimer
            flag[k] = false
        end
        k = k+1
    end

    cd(savePwd)

    # extrait les noms valides et retourne le vecteur correspondant
    finstances = allfiles[flag]
    return finstances
end



# ==============================================================================
# Parse an instance of 2SPA problem (bi-objective partionning problem) from a file

function parse2SPA(fname::String)

    f = open(fname)    
    nbctr, nbvar = parse.(Int, split(readline(f))) # nombre de contraintes , nombre de variables
    A = zeros(Int, nbctr, nbvar)                   # matrice des contraintes
    C = zeros(Int, 2,nbvar)                        # matrice des couts
    nb = zeros(Int, nbvar)
    for j in 1:nbvar
        flag = 1
        for valeur in split(readline(f))
            if flag == 1
                C[1,j] = parse(Int, valeur)
                flag +=1
            elseif flag == 2
                C[2,j] = parse(Int, valeur)
                flag +=1
            elseif flag == 3
                nb[j] = parse(Int, valeur)
                flag +=1
            else
                i = parse(Int, valeur)
                A[i,j] = 1
            end
        end
    end
    close(f)
    return C, A
end

