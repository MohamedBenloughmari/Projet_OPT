# =====================================================================
# 5OD21 - Optimisation Discrète - Projet "Course d'avions"
#
# Deux formulations PLNE inspirées du TSP orienté :
#   1) MTZ  : nombre polynomial de contraintes (Miller-Tucker-Zemlin)
#   2) CUT  : nombre exponentiel de contraintes de connexité,
#             ajoutées à la volée (lazy constraints) => Branch-and-Cut
#
# Utilisation :
#   julia course_avions.jl                          -> toutes les instances
#   julia course_avions.jl instances/aerodrome_20_1.txt
#   SOLVER=Gurobi julia course_avions.jl            -> avec Gurobi
# =====================================================================

using JuMP
const SOLVER = get(ENV, "SOLVER", "SCIP")
if SOLVER == "Gurobi"
    using Gurobi
    const GRB_ENV = Gurobi.Env()
    new_model() = Model(() -> Gurobi.Optimizer(GRB_ENV))
else
    using SCIP
    new_model() = Model(SCIP.Optimizer)
end

include("lecture_distances.jl")

const TIME_LIMIT = 600.0   # secondes par résolution
const EPS = 1e-6

# ---------------------------------------------------------------------
# Graphe : on garde uniquement les arcs (i,j) réalisables
#   - D[i,j] <= R (rayon d'action)
#   - pas d'arc entrant dans le départ d, pas d'arc sortant de l'arrivée f
# ---------------------------------------------------------------------
function build_arcs(n, d, f, R, D)
    return [(i, j) for i in 1:n for j in 1:n
            if i != j && D[i, j] <= R && j != d && i != f]
end

# ---------------------------------------------------------------------
# Partie commune aux deux formulations :
#   x[i,j] = 1 si l'avion vole directement de i à j
#   y[i]   = 1 si l'aérodrome i est visité
# ---------------------------------------------------------------------
function base_model(n, d, f, Amin, Nr, regions, D, A)
    m = new_model()
    set_silent(m)
    set_time_limit_sec(m, TIME_LIMIT)

    @variable(m, x[A], Bin)
    @variable(m, y[1:n], Bin)

    out_arcs = [[(i, j) for (i, j) in A if i == k] for k in 1:n]
    in_arcs  = [[(i, j) for (i, j) in A if j == k] for k in 1:n]

    @objective(m, Min, sum(D[i, j] * x[(i, j)] for (i, j) in A))

    # départ et arrivée
    @constraint(m, sum(x[a] for a in out_arcs[d]) == 1)
    @constraint(m, sum(x[a] for a in in_arcs[f]) == 1)
    @constraint(m, y[d] == 1)
    @constraint(m, y[f] == 1)

    # conservation du flot : on entre et on sort une fois de chaque aérodrome visité
    for k in 1:n
        k in (d, f) && continue
        @constraint(m, sum(x[a] for a in in_arcs[k]) == y[k])
        @constraint(m, sum(x[a] for a in out_arcs[k]) == y[k])
    end

    # nombre minimal d'aérodromes visités
    @constraint(m, sum(y) >= Amin)

    # chaque région visitée au moins une fois
    for k in 1:Nr
        @constraint(m, sum(y[i] for i in regions[k]) >= 1)
    end
    return m
end

# =====================================================================
# Formulation 1 : MTZ (polynomiale)
#   u[i] = position de i dans le chemin ; x[i,j] = 1  =>  u[j] >= u[i] + 1
# =====================================================================
function model_mtz(n, d, f, Amin, Nr, regions, D, A)
    m = base_model(n, d, f, Amin, Nr, regions, D, A)
    x = m[:x]
    @variable(m, 1 <= u[1:n] <= n)
    @constraint(m, u[d] == 1)
    for (i, j) in A
        @constraint(m, u[j] >= u[i] + 1 - n * (1 - x[(i, j)]))
    end
    return m
end

# =====================================================================
# Formulation 2 : coupes de connexité (exponentielle)
#   Pour tout S ⊆ V \ {d} et tout m ∈ S :
#       sum_{i ∉ S, j ∈ S} x[i,j] >= y[m]
#   ("si on visite m, au moins un arc doit entrer dans S")
# =====================================================================

# Ajoute (ou renvoie) les coupes pour un ensemble S
function cut_expr(m, A, S)
    x, y = m[:x], m[:y]
    inS = falses(size(y, 1))
    inS[S] .= true
    return sum(x[(i, j)] for (i, j) in A if !inS[i] && inS[j]; init = AffExpr(0.0))
end

# ---- Séparation entière : composantes non reliées au départ ---------
# xv[(i,j)] ∈ {0,1}. On parcourt le chemin depuis d ; les aérodromes
# visités qui ne sont pas atteints forment des sous-tours.
function separate_integer(n, d, A, xv)
    succ = zeros(Int, n)
    for (i, j) in A
        xv[(i, j)] > 0.5 && (succ[i] = j)
    end
    reached = falses(n)
    k = d
    while k != 0 && !reached[k]
        reached[k] = true
        k = succ[k]
    end
    subtours = Vector{Vector{Int}}()
    seen = copy(reached)
    for s in 1:n
        (seen[s] || succ[s] == 0) && continue
        cyc = Int[]
        k = s
        while !seen[k]
            seen[k] = true
            push!(cyc, k)
            k = succ[k]
        end
        push!(subtours, cyc)
    end
    return subtours
end

# ---- Séparation fractionnaire : coupe minimale d -> m ----------------
# Max-flot d'Edmonds-Karp sur la matrice des capacités cap = x*.
# Renvoie (valeur du flot, ensemble S côté puits de la coupe min).
function min_cut(cap::Matrix{Float64}, s::Int, t::Int)
    n = size(cap, 1)
    res = copy(cap)
    flow = 0.0
    while true
        parent = zeros(Int, n)
        parent[s] = s
        queue = [s]
        while !isempty(queue) && parent[t] == 0
            u = popfirst!(queue)
            for v in 1:n
                if parent[v] == 0 && res[u, v] > EPS
                    parent[v] = u
                    push!(queue, v)
                end
            end
        end
        parent[t] == 0 && break
        # capacité résiduelle du chemin augmentant
        b = Inf
        v = t
        while v != s
            b = min(b, res[parent[v], v]); v = parent[v]
        end
        v = t
        while v != s
            u = parent[v]
            res[u, v] -= b; res[v, u] += b; v = u
        end
        flow += b
    end
    # S = sommets non atteignables depuis s dans le graphe résiduel
    reach = falses(n); reach[s] = true
    queue = [s]
    while !isempty(queue)
        u = popfirst!(queue)
        for v in 1:n
            if !reach[v] && res[u, v] > EPS
                reach[v] = true; push!(queue, v)
            end
        end
    end
    return flow, findall(.!reach)
end

# Pour chaque m avec y*[m] > 0 : si mincut(d,m) < y*[m], la coupe est violée
function separate_fractional(n, d, A, xv, yv)
    cap = zeros(n, n)
    for (i, j) in A
        cap[i, j] = xv[(i, j)]
    end
    cuts = Tuple{Vector{Int},Int}[]
    covered = falses(n)          # évite de recalculer pour un m déjà dans un S trouvé
    for k in 1:n
        (k == d || yv[k] < EPS || covered[k]) && continue
        val, S = min_cut(cap, d, k)
        if val < yv[k] - 1e-4
            push!(cuts, (S, k))
            covered[S] .= true
        end
    end
    return cuts
end

function model_cut(n, d, f, Amin, Nr, regions, D, A)
    m = base_model(n, d, f, Amin, Nr, regions, D, A)
    x, y = m[:x], m[:y]
    m.ext[:ncuts] = 0
    SOLVER == "Gurobi" || return m

    # Gurobi : callback "lazy" appelé sur chaque solution entière candidate
    function lazy_cb(cb_data)
        xv = Dict(a => callback_value(cb_data, x[a]) for a in A)
        # on ne traite que les solutions entières
        any(v -> 0.01 < v < 0.99, values(xv)) && return
        for S in separate_integer(n, d, A, xv)
            e = cut_expr(m, A, S)
            for k in S
                con = @build_constraint(e >= y[k])
                MOI.submit(m, MOI.LazyConstraint(cb_data), con)
                m.ext[:ncuts] += 1
            end
        end
    end
    set_attribute(m, MOI.LazyConstraintCallback(), lazy_cb)
    return m
end

# Résolution de la formulation exponentielle.
#  - Gurobi : un seul Branch-and-Cut (coupes ajoutées par le callback)
#  - SCIP   : SCIP.jl ne gère pas les lazy callbacks de JuMP, on fait donc
#             de la génération de contraintes "à l'extérieur" : on résout,
#             on cherche les sous-tours, on ajoute les coupes, on recommence.
# Renvoie (temps total, nombre total de nœuds, nombre d'itérations)
function solve_cut!(m, n, d, A)
    if SOLVER == "Gurobi"
        optimize!(m)
        return solve_time(m), node_count(m), 1
    end
    x, y = m[:x], m[:y]
    t, nodes, it = 0.0, 0, 0
    while true
        set_time_limit_sec(m, max(1.0, TIME_LIMIT - t))
        optimize!(m)
        t += solve_time(m); nodes += node_count(m); it += 1
        primal_status(m) == FEASIBLE_POINT || return t, nodes, it
        termination_status(m) == OPTIMAL || return t, nodes, it
        xv = Dict(a => value(x[a]) for a in A)
        subtours = separate_integer(n, d, A, xv)
        isempty(subtours) && return t, nodes, it
        for S in subtours
            e = cut_expr(m, A, S)
            for k in S
                @constraint(m, e >= y[k])
                m.ext[:ncuts] += 1
            end
        end
    end
end

# ---------------------------------------------------------------------
# Bornes de relaxation continue (à la racine)
# ---------------------------------------------------------------------
function lp_bound_mtz(n, d, f, Amin, Nr, regions, D, A)
    m = model_mtz(n, d, f, Amin, Nr, regions, D, A)
    relax_integrality(m)
    optimize!(m)
    return is_solved_and_feasible(m) ? objective_value(m) : NaN
end

# Pour la formulation exponentielle : algorithme de plans coupants sur la
# relaxation continue, avec séparation exacte par coupe minimale.
function lp_bound_cut(n, d, f, Amin, Nr, regions, D, A)
    m = base_model(n, d, f, Amin, Nr, regions, D, A)
    relax_integrality(m)
    x, y = m[:x], m[:y]
    it = 0
    while true
        optimize!(m)
        is_solved_and_feasible(m) || return NaN, it
        xv = Dict(a => value(x[a]) for a in A)
        yv = value.(y)
        cuts = separate_fractional(n, d, A, xv, yv)
        isempty(cuts) && return objective_value(m), it
        for (S, k) in cuts
            @constraint(m, cut_expr(m, A, S) >= y[k])
        end
        it += 1
    end
end

# ---------------------------------------------------------------------
# Lecture de la solution : chemin d -> f
# ---------------------------------------------------------------------
function extract_path(n, d, f, A, x)
    succ = zeros(Int, n)
    for (i, j) in A
        value(x[(i, j)]) > 0.5 && (succ[i] = j)
    end
    path = [d]
    while path[end] != f && succ[path[end]] != 0 && length(path) <= n
        push!(path, succ[path[end]])
    end
    return path
end

function check_solution(path, d, f, Amin, Nr, regions, R, D)
    ok = path[1] == d && path[end] == f && allunique(path)
    ok &= length(path) >= Amin
    ok &= all(D[path[k], path[k+1]] <= R for k in 1:length(path)-1)
    ok &= all(any(i in path for i in regions[k]) for k in 1:Nr)
    return ok
end

function node_count(m)
    try
        return MOI.get(m, MOI.NodeCount())
    catch
        return -1
    end
end

# ---------------------------------------------------------------------
# Résolution d'une instance avec les deux formulations
# ---------------------------------------------------------------------
function solve_instance(file)
    n, d, f, Amin, Nr, R, regions, coords, D = readInstance(file)
    A = build_arcs(n, d, f, R, D)
    name = basename(file)
    println("\n=== $name : n=$n, d=$d, f=$f, Amin=$Amin, K=$Nr, R=$R, |A|=$(length(A))")

    results = []
    for form in (:MTZ, :CUT)
        if form == :MTZ
            lp = lp_bound_mtz(n, d, f, Amin, Nr, regions, D, A)
            m = model_mtz(n, d, f, Amin, Nr, regions, D, A)
            optimize!(m)
            t, nodes = solve_time(m), node_count(m)
        else
            lp, _ = lp_bound_cut(n, d, f, Amin, Nr, regions, D, A)
            m = model_cut(n, d, f, Amin, Nr, regions, D, A)
            t, nodes, _ = solve_cut!(m, n, d, A)
        end
        st = termination_status(m)
        if primal_status(m) == FEASIBLE_POINT
            obj = objective_value(m)
            bound = objective_bound(m)
            path = extract_path(n, d, f, A, m[:x])
            ok = check_solution(path, d, f, Amin, Nr, regions, R, D)
        else
            obj, bound, path, ok = NaN, NaN, Int[], false
        end
        ncuts = form == :CUT ? m.ext[:ncuts] : 0
        println("  [$form] statut=$st  obj=$obj  borneLP=$(round(lp, digits=2))  " *
                "temps=$(round(t, digits=2))s  noeuds=$nodes  coupes=$ncuts  valide=$ok")
        println("        chemin ($(length(path)) aérodromes) : ", join(path, " -> "))
        push!(results, (instance = name, form = form, status = st, obj = obj,
                        bound = bound, lp = lp, time = t, nodes = nodes,
                        ncuts = ncuts, nvisit = length(path), valid = ok))
    end
    return results
end

function write_table(results, out)
    open(out, "w") do io
        println(io, "| Instance | Formulation | Statut | Optimum | Borne LP racine | Gap LP (%) | Temps (s) | Nœuds | Coupes lazy | Aérodromes |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|")
        for r in results
            gap = isnan(r.obj) || isnan(r.lp) ? NaN : 100 * (r.obj - r.lp) / r.obj
            println(io, "| $(r.instance) | $(r.form) | $(r.status) | $(r.obj) | ",
                    "$(round(r.lp, digits=2)) | $(round(gap, digits=1)) | ",
                    "$(round(r.time, digits=2)) | $(r.nodes) | $(r.ncuts) | $(r.nvisit) |")
        end
    end
    println("\nTableau écrit dans $out")
end

# ---------------------------------------------------------------------
function main(args)
    dir = joinpath(@__DIR__, "instances")
    files = if isempty(args)
        # tri par taille : aerodrome_20_1.txt -> 20
        size_of(f) = (parse(Int, replace(split(f, "_")[2], ".txt" => "")), f)
        joinpath.(dir, sort(filter(endswith(".txt"), readdir(dir)), by = size_of))
    else
        args
    end
    all_results = []
    for file in files
        append!(all_results, solve_instance(file))
    end
    write_table(all_results, joinpath(@__DIR__, "resultats_$(SOLVER).md"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
