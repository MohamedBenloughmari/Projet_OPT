# 5OD21 – Projet « Course d'avions »

Fichiers :

| Fichier | Rôle |
|---|---|
| `lecture_distances.jl` | lecture d'une instance (fourni) |
| `course_avions.jl` | les deux formulations PLNE, les bornes LP, la résolution et le tableau de résultats |
| `resultats_SCIP.md` | tableau généré automatiquement par le script |
| `instances/` | instances fournies |

Lancer :

```bash
julia course_avions.jl                                   # toutes les instances
julia course_avions.jl instances/aerodrome_0.txt         # une instance
SOLVER=Gurobi julia course_avions.jl                     # avec Gurobi (si installé)
```

Par défaut le script utilise **SCIP** (gratuit, déjà installé). Avec `SOLVER=Gurobi`
il utilise Gurobi, comme demandé dans le complément du sujet.

---

## 1. Le problème

On a $n$ aérodromes $V = \{1,\dots,n\}$, un départ $d$, une arrivée $f$, $K$ régions,
un nombre minimal d'aérodromes $A_{min}$ et un rayon d'action $R$. On note
$V_k \subseteq V$ l'ensemble des aérodromes de la région $k$. On cherche un
**chemin élémentaire** de $d$ à $f$ qui :

- visite au moins $A_{min}$ aérodromes (départ et arrivée compris) ;
- passe par au moins un aérodrome de chaque région $k = 1,\dots,K$ ;
- n'utilise que des sauts de longueur $d_{ij} \le R$ ;
- minimise la distance totale.

C'est un **TSP orienté avec sélection des sommets** (proche du *Prize-Collecting
TSP* / *Orienteering*) : contrairement au TSP, on ne visite pas tous les sommets,
et on cherche un chemin $d \to f$ au lieu d'un circuit.

### Graphe utilisé

On ne crée une variable que pour les arcs utiles :

$$
A = \{ (i,j) \in V^2 \;:\; i \neq j,\; d_{ij} \le R,\; j \neq d,\; i \neq f \}
$$

- $d_{ij} > R$ : le vol direct est impossible, l'arc est supprimé (la contrainte de
  rayon est donc traitée « gratuitement », sans contrainte dans le modèle) ;
- pas d'arc entrant dans $d$ ni sortant de $f$ : le chemin commence en $d$ et finit en $f$.

On note $\delta^+(k) = \{(k,j) \in A\}$ les arcs sortant de $k$ et
$\delta^-(k) = \{(i,k) \in A\}$ les arcs entrant dans $k$.

## 2. Variables et contraintes communes

- $x_{ij} \in \{0,1\}$ pour $(i,j) \in A$ : 1 si l'avion vole directement de $i$ à $j$ ;
- $y_i \in \{0,1\}$ pour $i \in V$ : 1 si l'aérodrome $i$ est visité.

$$
\min \sum_{(i,j) \in A} d_{ij}\, x_{ij}
$$

$$
\begin{aligned}
(1)\quad & \sum_{(d,j) \in \delta^+(d)} x_{dj} = 1 && \text{on part de } d \\
(2)\quad & \sum_{(i,f) \in \delta^-(f)} x_{if} = 1 && \text{on arrive en } f \\
(3)\quad & y_d = y_f = 1 \\
(4)\quad & \sum_{(i,k) \in \delta^-(k)} x_{ik} = y_k && \forall k \in V \setminus \{d,f\} \quad \text{(on entre une fois dans } k \text{ s'il est visité)} \\
(5)\quad & \sum_{(k,j) \in \delta^+(k)} x_{kj} = y_k && \forall k \in V \setminus \{d,f\} \quad \text{(on en sort une fois)} \\
(6)\quad & \sum_{i \in V} y_i \ge A_{min} && \text{nombre minimal d'aérodromes} \\
(7)\quad & \sum_{i \in V_k} y_i \ge 1 && \forall k \in \{1,\dots,K\} \quad \text{(chaque région est visitée)} \\
& x_{ij} \in \{0,1\},\; y_i \in \{0,1\}
\end{aligned}
$$

Les contraintes (1)–(5) disent que les arcs choisis forment **un chemin de $d$ à $f$
plus, éventuellement, des circuits disjoints** (des « sous-tours »). Comme dans le
TSP, il faut interdire ces sous-tours : sinon le solveur « visite » des aérodromes
avec de petits circuits isolés (par exemple $a \to b \to a$) pour satisfaire (6) et
(7) à moindre coût. Les deux formulations ne diffèrent que sur cette partie.

## 3. Formulation 1 : MTZ (nombre polynomial de contraintes)

Inspirée de Miller–Tucker–Zemlin. On ajoute une variable continue $u_i \in [1, n]$,
la « position » de $i$ dans le chemin :

$$
\begin{aligned}
(8)\quad & u_d = 1 \\
(9)\quad & u_j \ge u_i + 1 - n\,(1 - x_{ij}) && \forall (i,j) \in A \\
& 1 \le u_i \le n && \forall i \in V
\end{aligned}
$$

Si $x_{ij} = 1$, alors $u_j \ge u_i + 1$ : les positions augmentent strictement le
long des arcs. Un circuit $i_1 \to i_2 \to \dots \to i_p \to i_1$ donnerait en
sommant $u_{i_1} \ge u_{i_1} + p$, impossible. Si $x_{ij} = 0$, la contrainte
devient $u_j \ge u_i + 1 - n$, toujours vraie puisque $u \in [1, n]$.

**Taille** : $|A| + 2n$ variables, $O(|A|)$ contraintes $\Rightarrow$ polynomial.

**Défaut** : la relaxation continue est faible. Avec $x_{ij}$ fractionnaire, le
terme $n(1 - x_{ij})$ est grand et (9) n'impose presque rien : la borne LP est
souvent loin de l'optimum, et le Branch-and-Bound explore beaucoup de nœuds.

## 4. Formulation 2 : coupes de connexité (nombre exponentiel de contraintes)

On remplace (8)–(9) par les contraintes de connexité :

$$
(10)\quad \sum_{\substack{(i,j) \in A \\ i \notin S,\; j \in S}} x_{ij} \;\ge\; y_m
\qquad \forall S \subseteq V \setminus \{d\},\; \forall m \in S
$$

« Si on visite $m$, au moins un arc doit entrer dans tout ensemble $S$ qui contient
$m$ mais pas $d$ ». Autrement dit, $m$ est atteignable depuis $d$ : il n'y a plus de
sous-tour. Un sous-tour $C$ violerait (10) avec $S = C$ (aucun arc n'entre dans
$C$, alors que $y_m = 1$).

C'est l'équivalent, pour notre problème, des contraintes d'élimination de sous-tours
du TSP $\left(\sum_{i \notin S,\, j \in S} x_{ij} \ge 1\right)$, avec $y_m$ à droite
parce qu'un aérodrome n'est pas forcément visité. Grâce à (4)–(5), on peut aussi
l'écrire sous la forme GSEC (*Generalized Subtour Elimination Constraints*) :

$$
\sum_{\substack{(i,j) \in A \\ i,j \in S}} x_{ij} \;\le\; \sum_{i \in S} y_i - y_m
\qquad \forall S \subseteq V \setminus \{d\},\; \forall m \in S
$$

**Taille** : $2^{n-1}$ ensembles $S$, multipliés par au plus $n$ choix de $m$
$\Rightarrow$ exponentiel. On ne peut pas tout écrire : on part du modèle (1)–(7)
et on ajoute **seulement les contraintes violées**, au fur et à mesure. Il faut
donc résoudre le *problème de séparation*.

### Problème de séparation

> Étant donné une solution $(x^*, y^*)$ (entière ou fractionnaire) qui vérifie
> (1)–(7), trouver $S \subseteq V \setminus \{d\}$ et $m \in S$ tels que
> $$\sum_{i \notin S,\, j \in S} x^*_{ij} < y^*_m,$$
> ou prouver qu'il n'en existe pas.

**Cas fractionnaire (séparation exacte en temps polynomial).** On voit $x^*_{ij}$
comme la capacité de l'arc $(i,j)$. Pour un $m$ fixé, la quantité

$$
\min_{S \,:\, m \in S,\; d \notin S} \;\sum_{i \notin S,\, j \in S} x^*_{ij}
$$

est exactement la **coupe minimale $d \to m$**, égale au **flot maximal de $d$ à
$m$** (théorème max-flot/min-cut). Donc :

```
pour chaque m ≠ d avec y*_m > 0 :
    calculer le flot max de d à m avec capacités x*
    si flot max < y*_m :  S = sommets non atteignables depuis d dans le graphe
                          résiduel  →  la contrainte (10) pour (S, m) est violée
```

Cela coûte au plus $n$ calculs de flot max, donc c'est polynomial. Dans le code :
`min_cut` (Edmonds–Karp) et `separate_fractional`.

**Cas entier (plus simple).** Si $x^*$ est entier, les arcs choisis forment un
chemin $d \to f$ plus des circuits. On suit le chemin depuis $d$ ; chaque circuit
$C$ non atteint donne une contrainte violée avec $S = C$ et $m \in C$ (on les
ajoute pour tous les $m \in C$). Dans le code : `separate_integer`.

### Mise en œuvre

- **Avec Gurobi** : la séparation entière est appelée dans un *lazy constraint
  callback* (`MOI.LazyConstraintCallback`) ; on obtient un vrai Branch-and-Cut,
  en une seule résolution.
- **Avec SCIP** : SCIP.jl ne gère pas les lazy callbacks de JuMP. On fait donc la
  boucle classique de génération de contraintes : résoudre le PLNE → chercher les
  sous-tours → ajouter les coupes → recommencer jusqu'à ce qu'il n'y ait plus de
  sous-tour. La solution obtenue est alors optimale pour le modèle complet (elle est
  optimale pour un modèle relâché et réalisable pour le modèle complet). Le temps et
  le nombre de nœuds indiqués sont la somme sur toutes les itérations.

## 5. Borne de la relaxation continue à la racine

Le complément demande la borne par relaxation continue à la racine pour chaque
formulation. On remplace $x_{ij} \in \{0,1\}$ et $y_i \in \{0,1\}$ par
$x_{ij} \in [0,1]$ et $y_i \in [0,1]$ :

- **MTZ** : on résout directement le PL obtenu (`relax_integrality`).
- **Coupes** : la relaxation continue contient *toutes* les contraintes (10). On la
  calcule par **plans coupants** : résoudre le PL, séparer avec le flot max
  (cas fractionnaire), ajouter les coupes violées, recommencer jusqu'à ce que la
  séparation ne trouve plus rien. On obtient exactement la borne LP du modèle
  exponentiel.

Le tableau donne aussi le *gap LP* :

$$
\text{gap} = \frac{z^* - z_{LP}}{z^*} \times 100
$$

où $z^*$ est la valeur optimale et $z_{LP}$ la borne de relaxation continue.

## 6. Validation

Sur l'instance d'exemple du sujet (`aerodrome_0.txt`, 6 aérodromes), les deux
formulations trouvent la valeur **12** avec le parcours **1 → 3 → 2 → 5**, ce qui
correspond à la solution donnée dans l'énoncé. Le script vérifie aussi chaque
solution de façon indépendante (`check_solution`) : chemin de $d$ à $f$ sans
répétition, au moins $A_{min}$ aérodromes, chaque saut $\le R$, chaque région visitée.

## 7. Résultats (SCIP, limite de 600 s par résolution)

RESULTS_PLACEHOLDER

## 8. Observations

OBSERVATIONS_PLACEHOLDER
