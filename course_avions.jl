# Projet 5OD21 - Course d'avions
# Formulation 1 : polynomiale (variables u)
# Formulation 2 : coupes, ajoutées par séparation

using JuMP, SCIP
include("lecture_distances.jl")

const TMAX = 600.0   # limite de temps par formulation (s)

# contraintes communes aux deux formulations
function modele_commun(n, d, f, Amin, Nr, regions, D, A)
    m = Model(SCIP.Optimizer)
    set_silent(m)
    set_time_limit_sec(m, TMAX)

    @variable(m, x[A], Bin)      # x[(i,j)] = 1 si on vole de i a j
    @variable(m, y[1:n], Bin)    # y[i] = 1 si on passe par i

    @objective(m, Min, sum(D[i, j] * x[(i, j)] for (i, j) in A))

    @constraint(m, sum(x[(i, j)] for (i, j) in A if i == d) == 1)
    @constraint(m, sum(x[(i, j)] for (i, j) in A if j == f) == 1)
    @constraint(m, y[d] == 1)
    @constraint(m, y[f] == 1)
    for k in 1:n
        if k != d && k != f
            @constraint(m, sum(x[(i, j)] for (i, j) in A if j == k) == y[k])
            @constraint(m, sum(x[(i, j)] for (i, j) in A if i == k) == y[k])
        end
    end
    @constraint(m, sum(y) >= Amin)
    for k in 1:Nr
        @constraint(m, sum(y[i] for i in regions[k]) >= 1)
    end
    return m
end

function modele_polynomial(n, d, f, Amin, Nr, regions, D, A)
    m = modele_commun(n, d, f, Amin, Nr, regions, D, A)
    x = m[:x]
    @variable(m, 1 <= u[1:n] <= n)
    @constraint(m, u[d] == 1)
    for (i, j) in A
        @constraint(m, u[j] >= u[i] + 1 - n * (1 - x[(i, j)]))
    end
    return m
end

# somme des x sur les arcs qui entrent dans S
function entrant(m, A, S)
    dansS = falses(length(m[:y]))
    dansS[S] .= true
    return sum(m[:x][(i, j)] for (i, j) in A if !dansS[i] && dansS[j]; init = AffExpr(0))
end

# solution entiere : on suit le chemin depuis d, ce qui reste forme des sous-tours
function sous_tours(n, d, A, xv)
    suiv = zeros(Int, n)
    for (i, j) in A
        if xv[(i, j)] > 0.5
            suiv[i] = j
        end
    end
    vu = falses(n)
    k = d
    while k != 0 && !vu[k]
        vu[k] = true
        k = suiv[k]
    end
    tours = Vector{Int}[]
    for s in 1:n
        if !vu[s] && suiv[s] != 0
            c = Int[]
            k = s
            while !vu[k]
                vu[k] = true
                push!(c, k)
                k = suiv[k]
            end
            push!(tours, c)
        end
    end
    return tours
end

# flot max de s vers t (chemins augmentants)
# renvoie la valeur du flot et les sommets du cote de t dans la coupe min
function flot_max(c, s, t)
    n = size(c, 1)
    c = copy(c)
    val = 0.0
    while true
        pred = zeros(Int, n)
        pred[s] = s
        file = [s]
        while !isempty(file)
            u = popfirst!(file)
            for v in 1:n
                if pred[v] == 0 && c[u, v] > 1e-6
                    pred[v] = u
                    push!(file, v)
                end
            end
        end
        if pred[t] == 0
            return val, findall(pred .== 0)
        end
        b = Inf
        v = t
        while v != s
            b = min(b, c[pred[v], v])
            v = pred[v]
        end
        v = t
        while v != s
            c[pred[v], v] -= b
            c[v, pred[v]] += b
            v = pred[v]
        end
        val += b
    end
end

# solution fractionnaire : pour chaque m, flot max de d vers m avec capacites x
function separation_fractionnaire(n, d, A, x, y)
    c = zeros(n, n)
    for (i, j) in A
        c[i, j] = value(x[(i, j)])
    end
    coupes = []
    for k in 1:n
        if k != d && value(y[k]) > 1e-6
            v, S = flot_max(c, d, k)
            if v < value(y[k]) - 1e-4
                push!(coupes, (S, k))
            end
        end
    end
    return coupes
end

function resoudre_coupes(n, d, f, Amin, Nr, regions, D, A)
    t0 = time()

    # 1) relaxation continue + coupes => borne a la racine
    r = modele_commun(n, d, f, Amin, Nr, regions, D, A)
    relax_integrality(r)
    lp = NaN
    toutes = []
    while true
        optimize!(r)
        if !is_solved_and_feasible(r)
            break
        end
        lp = objective_value(r)
        coupes = separation_fractionnaire(n, d, A, r[:x], r[:y])
        if isempty(coupes)
            break
        end
        for (S, k) in coupes
            @constraint(r, entrant(r, A, S) >= r[:y][k])
        end
        append!(toutes, coupes)
    end

    # 2) PLNE en gardant ces coupes, on ajoute les sous-tours tant qu'il y en a
    m = modele_commun(n, d, f, Amin, Nr, regions, D, A)
    x, y = m[:x], m[:y]
    for (S, k) in toutes
        @constraint(m, entrant(m, A, S) >= y[k])
    end
    noeuds = 0
    while true
        set_time_limit_sec(m, max(1.0, TMAX - (time() - t0)))
        optimize!(m)
        noeuds += MOI.get(m, MOI.NodeCount())
        if termination_status(m) != OPTIMAL || time() - t0 >= TMAX
            break
        end
        xv = Dict(a => value(x[a]) for a in A)
        tours = sous_tours(n, d, A, xv)
        if isempty(tours)
            break
        end
        for c in tours, k in c
            @constraint(m, entrant(m, A, c) >= y[k])
        end
    end
    return m, lp, time() - t0, noeuds
end

function chemin(n, d, f, A, x)
    suiv = zeros(Int, n)
    for (i, j) in A
        if value(x[(i, j)]) > 0.5
            suiv[i] = j
        end
    end
    p = [d]
    while p[end] != f && suiv[p[end]] != 0 && length(p) <= n
        push!(p, suiv[p[end]])
    end
    return p
end

function afficher(nom, m, lp, t, noeuds, n, d, f, A, Amin, Nr, regions)
    p = primal_status(m) == FEASIBLE_POINT ? chemin(n, d, f, A, m[:x]) : Int[]
    ok = !isempty(p) && p[end] == f && allunique(p) && length(p) >= Amin &&
         all(any(i in p for i in regions[k]) for k in 1:Nr)
    z = ok ? objective_value(m) : NaN
    binf = objective_bound(m)
    ecart = 100 * (z - lp) / z
    println("  $nom : ", termination_status(m), "  z = $z  borne inf = ", round(binf, digits = 2),
            "  borne LP = ", round(lp, digits = 2), "  ecart = ", round(ecart, digits = 1), "%",
            "  temps = ", round(t, digits = 2), "s  noeuds = $noeuds")
    if ok
        println("      chemin (", length(p), ") : ", join(p, " "))
    end
end

instances = ["aerodrome_0", "aerodrome_20_1", "aerodrome_20_2", "aerodrome_20_3",
             "aerodrome_20_4", "aerodrome_30_1", "aerodrome_40_1", "aerodrome_50_1",
             "aerodrome_70_1", "aerodrome_80_1", "aerodrome_90_1", "aerodrome_100_1",
             "aerodrome_100_2", "aerodrome_120_1", "aerodrome_150_1"]
fichiers = isempty(ARGS) ? ["instances/$nom.txt" for nom in instances] : ARGS

for fic in fichiers
    n, d, f, Amin, Nr, R, regions, coords, D = readInstance(fic)
    # on enleve les arcs trop longs, ceux qui arrivent en d et ceux qui partent de f
    A = [(i, j) for i in 1:n for j in 1:n if i != j && D[i, j] <= R && j != d && i != f]
    println("\n", fic, " : n = $n, |A| = ", length(A))

    r = modele_polynomial(n, d, f, Amin, Nr, regions, D, A)
    relax_integrality(r)
    optimize!(r)
    lp = objective_value(r)

    m = modele_polynomial(n, d, f, Amin, Nr, regions, D, A)
    optimize!(m)
    afficher("polynomiale", m, lp, solve_time(m), MOI.get(m, MOI.NodeCount()),
             n, d, f, A, Amin, Nr, regions)

    m, lp, t, noeuds = resoudre_coupes(n, d, f, Amin, Nr, regions, D, A)
    afficher("coupes", m, lp, t, noeuds, n, d, f, A, Amin, Nr, regions)
end
