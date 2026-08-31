include("model.jl")
using CSV
using DataFrames
using JLD2

### Input from shell skript ###
path = "./simulation_data/model_a/"
reps = 50
parameter = "q"
values = "range(0.1,step=0.1,stop=0.9)"
ID = parse(Int64,ARGS[5])
value = values[1] # this does a parameter scan just for the first q value,
# highly recommend to parallelize on server via slurm


### Parameter and initial values###

simulation_type_tup = read_sim_options(path*"simulation_options_new.txt")
par = CSV.File("new_parameter.csv") |> Dict

cell_ms = load("./initial_organoids/cell_ms_3d.jld2","cell_ms")
cell_ms.r0 = cell_ms.r
cell_ms.t0 = zeros(length(cell_ms.r))

connectivity = delaunay(cell_ms.xy')


### Second par change ###

dt = 0.01
par["maxerror"] = 1e-9

steady_state_parameter_scan(par,parameter,value,simulation_type_tup,cell_ms,
connectivity,path=path,replicates=reps,dt=dt,maxerror=par["maxerror"])