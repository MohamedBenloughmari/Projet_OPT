# ---------------------------------------------------------------
# Function that reads instances and returns formatted data
# Usage : n,d,f,Amin,Nr,R,regions,coords = readInstance("path/to/instance.txt")
# To use it in other files, use include("lecture_distances.jl")
# Parameter meanings :
# n : number of airports
# d : "depart"/start
# f : "fin"/end
# Amin : minimum number of airports to go through
# Nr : number of regions (all must be visited !)
# R : maximum range of a plane
# regions : array such that regions[region number] =
# list of airports in that region
# coords : coordinates of all airports
# ---------------------------------------------------------------

function readInstance(file::String)

  #Open file, extract data, close file
  o = open(file,"r")
  data = readlines(o)
  close(o)

  #Get simple numeric data
  n = parse(Int64,data[1])
  d = parse(Int64,data[2])
  f = parse(Int64,data[3])
  Amin = parse(Int64,data[4])
  Nr = parse(Int64,data[5])
  R = parse(Int64,data[9])

  #Get region definitions
  regions = Dict()
  for i in 1:Nr
    regions[i] = []
  end

  region_data = split(data[7]," ")
  for i in 1:n
    k = parse(Int64,region_data[i])
    if (k > Nr)
      # ex : aerodrome_90_1.txt déclare 6 régions mais utilise la région 7.
      # On respecte K : l'aérodrome est traité comme n'appartenant à aucune région.
      @warn "aérodrome $i : région $k > K=$Nr, ignorée"
    elseif (k != 0)
      append!(regions[k],i)
    end
  end

  #Get airport coordinates
  coords = Matrix{Int64}(zeros(n,2))
  for i in 1:n
    line = split(data[10+i]," ")
    coords[i,1] = parse(Int64,line[1])
    coords[i,2] = parse(Int64,line[2])
  end

  #Produce distance matrix
D = Matrix{Int64}(zeros(n,n))
for i in 1:n
  for j in 1:n
    D[i,j] = Int64(floor(sqrt((coords[j,1] - coords[i,1])^2
    + (coords[j,2] - coords[i,2])^2)))
  end
end

  return n,d,f,Amin,Nr,R,regions,coords,D
end