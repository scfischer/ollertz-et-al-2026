include("model.jl")
using CSV

path = "./simulation_data/model_a/"
reps = 50
parameter = "q"
values = "range(0.1,step=0.1,stop=0.9)"
ID = 1
value = values[ID] # this does a parameter scan just for the first q value,
# highly recommend to parallelize on server via slurm

### Paramater ###
simulation_type_tup = read_sim_options("simulation_options_new.txt")
par = CSV.File("new_parameter.csv") |> Dict

### Time ###
time_ss = Time_struct(5000,1000)


### Initial conditions
cells_dict = Dict(
    "xy" => ([-0.1 0 0.3; 0.1 0 0.3; 0.05 0.05 0.3; 0.2 0.2 0.2; -0.2 -0.2 -0.2]),
    "r" => [0.9, 0.9, 0.9, 0.9, 0.9],
    "u" => [0.75, 0.75, 0.75, 0.75, 0.75],
    "v" => [0.75, 0.75, 0.75, 0.75, 0.75],
)
cell_ms = Cell_struct(cells_dict)

### time course data ###
names = (:xy, :r, :u, :v)
data_tup = create_data_tup(names)

parameter_scan(par,parameter,value,simulation_type_tup,cell_ms,time_ss,data_tup,path=path,replicates=reps,three_d=true)