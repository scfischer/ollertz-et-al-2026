using Random, Distributions #used in cell_division
using LinearAlgebra # used in Forces, non_random_cell_devision,signalling_distance
using Distances # used in Forces
using Serialization # used in save_data
using Combinatorics
using Graphs
using Plots,GraphRecipes
using MiniQhull
using ColorSchemes 
using Measures



function read_sim_options(file::String)
    return eval(Meta.parse(read(file,String)))
end

struct Time_struct
    T::Int64
    N::Int64
    dt::Float64
    t::Array{Float64}
end

function Time_struct(T::Int64,N::Int64)

    t = collect(range(0,T,length=N*T))   # Time at each timestep
    dt = t[2]-t[1]           # timestep
    return Time_struct(T,N*T,dt,t)
end

mutable struct Cell_struct
    xy::Array{Float64}
    r::Array{Float64}
    u::Array{Float64}
    v::Array{Float64}
    r0::Array{Float64}
    t0::Array{Float64}
    s::Array{Float64}
end

function Cell_struct(xy::Array{Float64},r::Array{Float64},u::Array{Float64},v::Array{Float64})

    r0 = r
    t0 = zeros(length(r))
    s = zeros(length(r))
    return Cell_struct(xy,r,u,v,r0,t0,s)
end

function Cell_struct(cell_dict::Dict{String, Array{Float64}})

    r0 = cell_dict["r"]
    t0 = zeros(length(cell_dict["r"]))
    s = zeros(length(cell_dict["r"]))
    return Cell_struct(cell_dict["xy"],cell_dict["r"],cell_dict["u"],cell_dict["v"],r0,t0,s)
end


function create_data_tup(names::NTuple{}=(:xy, :r, :u, :v,:s))

    values = []
    for i in 1:length(names)
        push!(values,Array{Array{Float64}}(undef,0))
    end
    data_tuple  = NamedTuple{names}(values)
    return  data_tuple
end

function create_data_tup(names::Symbol)

    data_tuple  = (names = Array{Array{Float64}}(undef,0),)
    return  data_tuple
end

function push_cells_to_data!(data::NamedTuple,cells::Cell_struct)

    for i in fieldnames(typeof(data))
        push!(getproperty(data,i),copy(getfield(cells,i))) # maybe copy is needed maybe not, jsut to be sure

    end
end


# needs to be changed to: r_t(r_max,k,t0,r0,t) = r_max./(1 .+ ((r_max - r0)/(r0)) .* exp.(-k * r_max * t))
function cell_growth(t, cell_ms, par)

    # Packages: none

    cell_ms.r = par["r_max"] .- exp.(-par["k"] .* (t .- cell_ms.t0)) .*(par["r_max"] .- cell_ms.r0)

    return cell_ms
end


function division_probability(r, r_max)

    # Packages: none

    c = 100 / r_max
    y = 0.95 * r_max
    r_min = 0.9 * r_max

    b = (1 + exp(c * (r_max - y))) * (1 + exp(c * (r_min - y)))/(exp(c * (r_min - y)) - exp(c * (r_max - y)))
    a = -b/(1 + exp(c*(r_min - y)))
    P = a+b/(1 + exp(c*(r - y)))

    return maximum([0, P])
end

function cell_division(cell_ms, r_old, t,dt, par,simulation_type_tup)

    # Packages: Random, Distributions


    # Choose where division will randomly occur    
    p_old = division_probability.(r_old, par["r_max"])
    p = division_probability.(cell_ms.r, par["r_max"])

    prob = (p .- p_old) ./ (1 .- p_old)
    prob[r_old .== cell_ms.r0]  = p[r_old .== cell_ms.r0]

        ####### test für dt ######
    prob = division_probability.(cell_ms.r, par["r_max"]) * dt

    random_numbers = rand(length(cell_ms.r))
    indices = findall(x -> x < 0, random_numbers - par["cell_div_scaling"]*prob)


    if !isempty(indices)

         # New radius based on the area of the mother cell being two times that of the daughter cells
        r_new = cell_ms.r[indices]/2^(1/2)

        # distance between the two daughter cells
        distance = rand.(Normal.((cell_ms.r[indices] - r_new)/2,0.1*(cell_ms.r[indices] - r_new)/2))

        # kind of cell devision
        cell_ms.xy = simulation_type_tup.cell_div_func(cell_ms.xy, indices, r_new, distance,simulation_type_tup)



        # radii update of mothers cells and adding radii of daughter cells # do the radii really get updated this way?
        for i in range(1,length(r_new),step=1)
            push!(cell_ms.r, r_new[i])
            push!(cell_ms.r0, r_new[i])
        end
        
        # t0 update of mothers cells and adding t0 of daughter cells
        # t0 is used to calculate cell growth
        cell_ms.t0[indices] .= t
        for i in range(1,length(cell_ms.t0[indices]),step=1)
            push!(cell_ms.t0, cell_ms.t0[indices[i]])
        end
        
        if !isempty(cell_ms.u)
            # G update of mothers cells and adding t0 of daughter cells
            # G is used to calculate cell growth
            if simulation_type_tup.tf_split
                cell_ms.v[indices] .= cell_ms.v[indices]/2
                for i in range(1,length(cell_ms.v[indices]),step=1)
                    push!(cell_ms.v, cell_ms.v[indices[i]]/2)
                end

                # N update of mothers cells and adding t0 of daughter cells
                # N is used to calculate cell growth
                cell_ms.u[indices] .=cell_ms.u[indices]/2
                for i in range(1,length(cell_ms.u[indices]),step=1)
                    push!(cell_ms.u, cell_ms.u[indices[i]]/2)
                end
            else
                cell_ms.v[indices] .= cell_ms.v[indices]
                for i in range(1,length(cell_ms.v[indices]),step=1)
                    push!(cell_ms.v, cell_ms.v[indices[i]])
                end

                # N update of mothers cells and adding t0 of daughter cells
                # N is used to calculate cell growth
                cell_ms.u[indices] .=cell_ms.u[indices]
                for i in range(1,length(cell_ms.u[indices]),step=1)
                    push!(cell_ms.u, cell_ms.u[indices[i]])
                end
            end
        end
    end

    return cell_ms
end

function random_angle_cell_devision(xy,indices,r_new,distance,simulation_type_tup)

    # Packages: Random

    # angle_value of cell division
    angle_value = rand(length(r_new))*2*π

    # displacement vectors
    dx = distance .* cos.(angle_value)
    dy = distance .* sin.(angle_value)

    # displacement
    dxy = [dx;dy]
    dxy = reshape(dxy,(2,Int(length(dxy)/2))) 
    xy1 = xy[indices,:] + dxy'
    xy2 = xy[indices,:] - dxy'

    # position update of mother cells
    xy[indices,:] = xy1

    # adding positions of daughter cells
    xy = vcat(xy,xy2)

    return xy

end

function non_random_cell_devision(xy,indices,r_new, distance,simulation_type_tup)

    # Packages: LinearAlgebra

    # get center of organoid
    center_xy = round.(mean(eachrow(xy)),digits= 3)
   
    # define distance vector to center for each mother cell
    dxy = xy[indices,:]  .- center_xy'  

    # scale it to length 1
    unit_dxy = dxy ./sqrt.(sum(eachcol(dxy.^2)))

    if simulation_type_tup.func_orth == true # if false, cell divide along axis from mother cell to organoid center, if true orthogonal to this axis 
        for i in range(1,length(dxy[:,1])) # nullspace. does not work -.-
            dxy[i,:] = nullspace(dxy[i,:]')
        end
    end

    # displacement
    xy1 = xy[indices,:] .+ distance .* unit_dxy
    xy2 = xy[indices,:] .- distance .* unit_dxy
    
    # position update of mother cells
    xy[indices,:] = xy1

    # adding positions of daughter cells
    xy = vcat(xy,xy2)

    return xy

end
function forces(cell_ms,xy, r, par,simulation_type_tup)

    #Packages: Distances, LinearAlgebra

    # alpha is the stiffness
    # sigma determines where the cells are ideally positioned in terms of their radii
    # F0 is a scaling factor

    distance = pairwise(Euclidean(),xy, dims = 1)

    r_pairwise = r' .+ r
    x_pairwise = xy[:,1]' .- xy[:,1]
    y_pairwise = xy[:,2]' .- xy[:,2]
    
    
    if simulation_type_tup.affinity
        # calculate affinity
        eta_a = ones(length(cell_ms.v),length(cell_ms.v))
        for k in range(1,length(cell_ms.r))
            for i in range(1,length(cell_ms.r))
                eta_a[k,i] = ((cell_ms.u[k] + cell_ms.v[k]) * (cell_ms.u[i] + cell_ms.v[i]))/(cell_ms.v[k] + cell_ms.v[i] + par["beta"]) + 1/par["beta"]
            end
        end
        # Absolute values of forces according to Morse potential
        F = par["F0"]*2*par["alpha"]*(exp.(-2*par["alpha"]*(distance-r_pairwise*par["sigma"])) - eta_a .* exp.(-par["alpha"]*(distance-r_pairwise*par["sigma"])))

    else

        # Absolute values of forces according to Morse potential
        F = par["F0"]*2*par["alpha"]*(exp.(-2*par["alpha"]*(distance-r_pairwise*par["sigma"])) - exp.(-par["alpha"]*(distance-r_pairwise*par["sigma"])))
    end
    #F[dist > r_pairwise] = 0

    # Fill distance matrix with inf on diagonal
    distance[diagind(distance)].=Inf32

    # x- and y-direction of forces
    Fx = F.*(x_pairwise)./distance
    Fy = F.*(y_pairwise)./distance

    # Sum of all forces acting on each cell as a vector
    force = hcat(sum(Fx, dims = 2),sum(Fy, dims = 2))
return (force)

end


function displacement(cell_ms, dt, par,simulation_type_tup)

    # Packages: none

    F = forces(cell_ms, cell_ms.xy, cell_ms.r, par,simulation_type_tup)

    if simulation_type_tup.brownian
        cell_ms.xy = cell_ms.xy - dt*F  + par["tD"] * rand(Normal(0, dt^(1/2)), length(cell_ms.xy[:,1]),2)
    else
        cell_ms.xy = cell_ms.xy - dt*F
    end

    return cell_ms
end

function graph_from_delaunay(connectivity,xy)

    # Packages: Graphs, Combinatorics

    g = SimpleGraph(zeros(length(xy[1,:]),length(xy[1,:]))) # define graph with number of notes as cells and no edges

    for k in 1:length(connectivity[1,:])    # iterate over all triangles from delaunay tesselation
        # get all edges of one triangle:
        edge_list = combinations(connectivity[:,k],2) |>collect # same as collect(combinations(connectivity[:,k],2))
        
        for i in 1:length(edge_list)    # iterate over all edges
            add_edge!(g,edge_list[i][1],edge_list[i][2])    # add edges to graph
        end
    end
return(g)
end

function matrix_from_delaunay(connectivity,cell_ms,simulation_type_tup)
    xy = cell_ms.xy'
    # Packages: Combinatorics

    m = zeros(length(xy[1,:]),length(xy[1,:])) # define graph with number of notes as cells and no edges
    for k in 1:length(connectivity[1,:])    # iterate over all triangles from delaunay tesselation
        # get all edges of one triangle:
        edge_list = combinations(connectivity[:,k],2) |>collect # same as collect(combinations(connectivity[:,k],2))

        if simulation_type_tup.check_distance # check if distance between to points is greater than the sum of their radii
            for i in 1:length(edge_list)    # iterate over all edges
                if euclidean(xy[:,edge_list[i][1]],xy[:,edge_list[i][2]]) < (cell_ms.r[edge_list[i][1]] + cell_ms.r[edge_list[i][2]])
                    m[edge_list[i][1],edge_list[i][2]] = 1    # add edges to graph
                    m[edge_list[i][2],edge_list[i][1]] = 1 
                end
            end
        else
            for i in 1:length(edge_list)    # iterate over all edges
                m[edge_list[i][1],edge_list[i][2]] = 1    # add edges to graph
                m[edge_list[i][2],edge_list[i][1]] = 1 
            end
        end
    end
return(m)
end

function dijkstra_from_matrix(m,source) #used in signalling_distance

    # Packages: none
    
    not_visited = fill(true,1,length(m[1,:]));
    dist = fill(1000,1,length(m[1,:]));
    connected = findall(x-> x > 0,m[source,:]);
    dist[connected] .= 1
    dist_copy = deepcopy(dist)
    dist[source] = 0

    # 1 step

    not_visited[source] = false
    counter = 0
    while sum(not_visited) > 0 && counter <1000

        i = findmin(dist_copy)[2][2]

        if not_visited[i] ==true  
            connected_new = findall(x-> x > 0,m[i,:])
    
            for j in connected_new
                if   dist[j] > 1 + dist[i] || dist[j] == 1000
                    dist[j] = 1 + dist[i]
                    dist_copy[j] = dist[j] + dist[i]
       
                end
            end
            not_visited[i] = false
        end 
        dist_copy[i] = 1000
        counter +=1
    end
    
    return (dist)
end


function remove_diagonal(x) #used in signalling_distance

        # Packages: none
    
    matrix = Array{Float64}(undef, size(x, 1), size(x, 1) - 1)
    for i = 1:size(matrix, 1)
        for j ∈ 1:size(matrix, 2)
            if i > j
                matrix[i, j] = x[i, j]
            else
                matrix[i, j] = x[i, j+1]
            end
        end    
    end
    return matrix
end

function add_diagonal(x) #used in signalling_distance

        # Packages: none
    
    matrix = zeros(size(x, 1), size(x, 1))
    for i = 1:size(matrix, 1)
        for j ∈ 1:size(x, 2)
            if i > j
                matrix[i, j] = x[i, j]
            else
                matrix[i, j+1] = x[i, j]
            end
        end    
    end
    return matrix
end


function signalling_distance(cell_ms,adj_matrix,par,simulation_type_tup) # distance based signalling

    # Packages: LinearAlgebra
    
   infl_matrix = calculate_infl_matrix(adj_matrix,par,simulation_type_tup)

   s = zeros(length(cell_ms.u))
   mul!(s,infl_matrix,cell_ms.u)
   
   cell_ms.s = s
   return(cell_ms)
end

function calculate_infl_matrix(adj_matrix,par,simulation_type_tup) # distance based signalling

    # Packages: LinearAlgebra
    
   distance = zeros(dim(adj_matrix),dim(adj_matrix))
   for i in range(1,dim(adj_matrix))
       distance[i,:]= dijkstra_from_matrix(adj_matrix,i)
   end

   distance = remove_diagonal(distance)
   infl_matrix = zeros(length(distance[:,1]),length(distance[1,:]))

   # Scaling Sascha (???) 

   """
   for  i in range(1,length(distance[:,1]))
       infl_matrix[i,:] = par["q"] .^ (distance[i,:].-1) * (maximum(sum(par["q"] .^(distance[i,:].-1))))^(-1)
   end
    """

    # Scaling Simon (Topography):
   if simulation_type_tup.topography 
        scaling = (maximum(sum(par["q"] .^(distance.-1),dims=2)))^(-1)
        for  i in range(1,length(distance[:,1]))
            infl_matrix[i,:] = par["q"] .^ (distance[i,:].-1) * scaling
        end
    
    else 
        
        for  i in range(1,length(distance[:,1]))
            infl_matrix[i,:] = par["q"] .^ (distance[i,:].-1) * (maximum(sum(par["q"] .^(distance[i,:].-1))))^(-1)
        end
    end 


   infl_matrix = round.(add_diagonal(infl_matrix),digits = 4)

    if simulation_type_tup.donator 
        infl_matrix = infl_matrix'
    end

   return infl_matrix
    
end


function signalling_nn(cell_ms, adj_matrix,par,simulation_type_tup)  #nearest neighbour signalling

    # Packages: none

    cell_ms.s = Array{Float64}(undef,0)
    for i in range(1,length(cell_ms.u))
        push!(cell_ms.s,dot(cell_ms.u,adj_matrix[i,:]))
    end
    return cell_ms
end


function transcription(cell_ms,dt,par,simulation_type_tup) 

    # Packages: none
   
    a = exp(-par["eps_G"])
    b = exp(-par["eps_N"])
    c = exp(-par["eps_S"])
    d = exp(-par["eps_GS"])


    pu =        (b*cell_ms.u)       ./ (1 .+ a*cell_ms.v .* (1 .+ d*c*cell_ms.s) + b*cell_ms.u + c*cell_ms.s)
    pv = (a*cell_ms.v) .* (1 .+ d*c*cell_ms.s) ./ (1 .+ a*cell_ms.v .* (1 .+ d*c*cell_ms.s) + b*cell_ms.u + c*cell_ms.s)
                        
    cell_ms.u = cell_ms.u .+ par["tau"] * dt *par["dt_trans"] .* (par["r_N"]*pu - par["Gamma_N"]*cell_ms.u)
    cell_ms.v = cell_ms.v .+ par["tau"] * dt *par["dt_trans"] .* (par["r_G"]*pv - par["Gamma_G"]*cell_ms.v)

    if simulation_type_tup.set_to_zero 
        cell_ms.u[cell_ms.u .< par["trans_tol"]] .= 0
        cell_ms.v[cell_ms.v .< par["trans_tol"]] .= 0
    end

    return cell_ms
end

### Data handling, plotting###

function transform_data_to_matrix(data)

    # Packages: none

    results_matrix = Array{Union{Missing, Float64}}(missing,(length(data),length(data[end])))
    for j in range(1,length(data))
        for i in range(1,length(data[j]))
            results_matrix[j,i] = data[j][i]
        end 
    end
    return results_matrix
end

function DrawTri!(T) """ from Alphashapes.jl """

    # Packages: non, but needs input from AlphaShapes - replaced by: draw_edges_from_matrix!(), draw_edges_from_connectivity!()

    plot!([T[1,1],T[2,1]],[T[1,2],T[2,2]],label="",alpha= 0.4,color="gray")
    plot!([T[2,1],T[3,1]],[T[2,2],T[3,2]],label="",alpha= 0.4,color="gray")
    plot!([T[3,1],T[1,1]],[T[3,2],T[1,2]],label="",alpha= 0.4,color="gray")
end

function draw_edges_from_connectivity!(connectivity,xy)

    # Packages: Combinatorics, Plots

    for k in 1:length(connectivity[1,:])    # iterate over all triangles from delaunay tesselation
        # get all edges of one triangle:
        edge_list = combinations(connectivity[:,k],2) |>collect # same as collect(combinations(connectivity[:,k],2))
        
        for i in 1:length(edge_list)    # iterate over all edges

            a = [xy[edge_list[i][1],1],xy[edge_list[i][2],1]]
            b = [xy[edge_list[i][1],2],xy[edge_list[i][2],2]]
            
            plot!(a,b,label="",alpha= 0.4,color="gray") #plot edges
        end
    end
end

function draw_edges_from_matrix!(matrix,xy)

    # Packages: Plots

    # matrix = connectivity_matrix xy = cell positions

    for k in range(1,length(matrix[1,:])) # iterate over matrix
        for i in range(1,k-1)
            if matrix[k,i] == 1 # if cells i and k are connected, a line is drawn between their positions
                a = [xy[i,1],xy[k,1]]
                b = [xy[i,2],xy[k,2]]
                plot!(a,b,label="",alpha= 0.4,color="gray")
            end
        end
    end
end

function circleShape(h,k,r)

    # Packages: none

    θ = LinRange(0,2*π,500)
    h .+ r*sin.(θ), k .+ r*cos.(θ)
end

function draw_radius!(cell_ms ,lw = 0.5, c = :blue, linecolor = :black, fillalpha = 0.1)

    # Packages: Plots

    plot!(circleShape.(cell_ms.xy[:,1],cell_ms.xy[:,2],cell_ms.r),seriestype = [:shape], lw = lw, c = c, linecolor = linecolor ,fillalpha=fillalpha,label="",aspect_ratio=1)

end


### save simulation details and data ###

function save_data(data,name,path="")
    serialize(path * name * ".jld2",data)
    
end

function save_sim_type(sim_type,name,path="",type=".txt",overwrite=false)
    open(path*name*type,"a") do io 
        if filesize(path*name*type) == 0 
            println(io,"simulation_type_tup = ",sim_type)
        else 
            println("file already exist")
            if overwrite
                println(io,"simulation_type_tup = ",sim_type)
            end
        end
    end
end

function save_par(par_dict,name,path="")
    CSV.write(path*name*".csv",par_dict)
end

function run_sim(cell_ms,time_ss,par,data_tup,simulation_type_tup)

    if simulation_type_tup.anim
        a = Animation()
    end

    if simulation_type_tup.transcr || simulation_type_tup.draw_edges
        connectivity = delaunay(cell_ms.xy')
        matrix = matrix_from_delaunay(connectivity, cell_ms, simulation_type_tup)
    end


    for i in range(1,time_ss.N,step=1)

        if par["max_cells"] <= length(cell_ms.r) && simulation_type_tup.cell_number_check
            break
        end

        ### physics simulation ###
        if simulation_type_tup.physics
            r_old = cell_ms.r
            cell_ms = cell_growth(time_ss.t[i],cell_ms, par)
            cell_ms = displacement(cell_ms, time_ss.dt, par,simulation_type_tup)
        end

        ### cell division ###
        if simulation_type_tup.cell_division
            cell_ms = cell_division(cell_ms, r_old, time_ss.t[i],time_ss.dt, par,simulation_type_tup)
        end

        ### graph ###
        if simulation_type_tup.transcr || simulation_type_tup.draw_edges
            connectivity = delaunay(cell_ms.xy')
            matrix = matrix_from_delaunay(connectivity, cell_ms, simulation_type_tup)
        end

        ### transcription ###
        if simulation_type_tup.transcr    
            cell_ms  = simulation_type_tup.sign_func(cell_ms, matrix,par,simulation_type_tup)    
            cell_ms  = transcription(cell_ms,time_ss.dt,par,simulation_type_tup) 

        end

        ### keep timecourse data ###
        if simulation_type_tup.save_data && i % par["save_intervall"] == 0
            
            push_cells_to_data!(data_tup,cell_ms)

        end
   
        ### plotting ###
        if simulation_type_tup.anim
            p = scatter(cell_ms.xy[:,1],cell_ms.xy[:,2],label="",xlabel="x",ylabel="y",aspect_ratio=:equal,marker_z=cell_ms.u,clim = (0, 0.1),title = "Organoid at time: " * string(round(time_ss.t[i],digits=1)),
                colorbar_title="u")
            if simulation_type_tup.draw_edges
                draw_edges_from_matrix!(matrix,cell_ms.xy)
            end
            if simulation_type_tup.draw_radius
                draw_radius!(cell_ms)
            end
            frame(a,p)
        end
        if simulation_type_tup.track_progress
            if i in Int.(time_ss.N .* [0.2,0.4,0.6,0.8])         
                println("##")
                
            end
        end
        
    end

    results =[]
    push!(results,cell_ms)

    if simulation_type_tup.save_data
            
        push!(results,data_tup)
    end
    
    if simulation_type_tup.anim
        push!(results,a)

    end
    return(results)
   
end

function run_steady_state(cell_ms, matrix, par,simulation_type_tup, dt, maxsteps, maxerror=1e-6;save_time_course=false)

    # calculate the influence matrix
    infl_matrix = calculate_infl_matrix(matrix, par,simulation_type_tup)
    
    # define empty array of floats    
    error = zeros(Float64,maxsteps)  

    if save_time_course
        names = (:xy, :r, :u, :v, :s)
        data_tup = create_data_tup(names)
        push_cells_to_data!(data_tup,cell_ms)
    end
    
    for i in range(1,maxsteps)
    
        # copy the current state of the system
        control_u = copy(cell_ms.u)
        control_v = copy(cell_ms.v)
    
        # calculate the signalling 
        #cell_ms  = simulation_type_tup.sign_func(cell_ms,infl_matrix,par,simulation_type_tup)    
          
        s = zeros(length(cell_ms.u))
        mul!(s,infl_matrix,cell_ms.u) 
        cell_ms.s = s
        
        # calculate the transcription
        cell_ms  = transcription(cell_ms, dt, par,simulation_type_tup)
        if save_time_course
            push_cells_to_data!(data_tup,cell_ms)
        end
        # calculate the error
    
        error[i] = norm(cell_ms.u - control_u) + norm(cell_ms.v - control_v)
        if error[i] < maxerror
            println(i)
            break
        end
    end 
        
        if save_time_course && i % par["save_intervall"] == 0
            return(cell_ms,data_tup)
        else
            return(cell_ms)
        end
    end
### 3D ###

function run_sim_3d(cell_ms,time_ss,par,data_tup,simulation_type_tup)

    if simulation_type_tup.anim
        a = Animation()
    end

    connectivity = delaunay(cell_ms.xy')
    matrix =[]

    for i in range(1,time_ss.N,step=1)
        
        if par["max_cells"] <= length(cell_ms.r) && simulation_type_tup.cell_number_check
            break
        end

        ### physics simulation ###
        if simulation_type_tup.physics
            r_old = cell_ms.r
            cell_ms = cell_growth(time_ss.t[i],cell_ms, par)
            cell_ms = displacement_3d(cell_ms, time_ss.dt, par,simulation_type_tup)
        end

        ### cell division ###
        if simulation_type_tup.cell_division
        cell_ms = cell_division_3d(cell_ms, r_old, time_ss.t[i],time_ss.dt, par,simulation_type_tup)
        end

        ### graph ###
        if simulation_type_tup.transcr || simulation_type_tup.draw_edges
            connectivity = delaunay(cell_ms.xy')
            matrix = matrix_from_delaunay(connectivity, cell_ms, simulation_type_tup)

        end

        ### transcription ###
        if simulation_type_tup.transcr    
            cell_ms  = simulation_type_tup.sign_func(cell_ms, matrix, par,simulation_type_tup)
            cell_ms  = transcription(cell_ms,time_ss.dt,par,simulation_type_tup) 

        end

        ### keep timecourse data ###
        if simulation_type_tup.save_data  && i % par["save_intervall"] == 0

            push_cells_to_data!(data_tup,cell_ms)

        end
   
        ### plotting ###
        if simulation_type_tup.anim

            p = scatter(cell_ms.xy[:,1],cell_ms.xy[:,2],cell_ms.xy[:,3],label="",xlabel="x",ylabel="y",aspect_ratio=:equal,marker_z=cell_ms.u,clim = (0, 0.1),title = "Organoid at time: " * string(round(time_ss.t[i],digits=1)),
                colorbar_title="u")
            if simulation_type_tup.draw_edges
                draw_edges_from_connectivity_3d!(connectivity,cell_ms.xy)
            end
            frame(a,p)
        end

        ### track progress ###
        if simulation_type_tup.track_progress
     
            if i in Int.(time_ss.N .* [0.2,0.4,0.6,0.8])    

                println("##")
                
            end
        end
        
    end

    results =[]
    push!(results,cell_ms)

    if simulation_type_tup.save_data
            
        push!(results,data_tup)
    end
    
    if simulation_type_tup.anim
        push!(results,a)

    end
    return(results)
   
end



function displacement_3d(cell_ms, dt, par,simulation_type_tup)

    # Packages: none

    F = forces_3d(cell_ms,cell_ms.xy, cell_ms.r, par,simulation_type_tup)

    if simulation_type_tup.brownian

        cell_ms.xy = cell_ms.xy - dt*F  + par["tD"] * rand(Normal(0, dt^(1/2)), length(cell_ms.xy[:,1]),3)
    else

        cell_ms.xy = cell_ms.xy - dt*F
    end
  

    return cell_ms
end


function forces_3d(cell_ms,xy, r, par,simulation_type_tup)
    # alpha is the stiffness
    # sigma determines where the cells are ideally positioned in terms of their radii
    # F0 is a scaling factor
    distance = pairwise(Euclidean(),xy, dims = 1)

    r_pairwise = r' .+ r
    x_pairwise = xy[:,1]' .- xy[:,1]
    y_pairwise = xy[:,2]' .- xy[:,2]
    z_pairwise = xy[:,3]' .- xy[:,3]
    
    if simulation_type_tup.affinity
        # calculate affinity
        eta_a = ones(length(cell_ms.v),length(cell_ms.v))
        for k in range(1,length(cell_ms.r))
            for i in range(1,length(cell_ms.r))
                eta_a[k,i] = ((cell_ms.u[k] + cell_ms.v[k]) * (cell_ms.u[i] + cell_ms.v[i]))/(cell_ms.v[k] + cell_ms.v[i] + par["beta"]) + 1/par["beta"]
            end
        end
        # Absolute values of forces according to Morse potential
        F = par["F0"]*2*par["alpha"]*(exp.(-2*par["alpha"]*(distance-r_pairwise*par["sigma"])) - eta_a .* exp.(-par["alpha"]*(distance-r_pairwise*par["sigma"])))

    else

        # Absolute values of forces according to Morse potential
        F = par["F0"]*2*par["alpha"]*(exp.(-2*par["alpha"]*(distance-r_pairwise*par["sigma"])) - exp.(-par["alpha"]*(distance-r_pairwise*par["sigma"])))
    end
    # Fill distance matrix with inf on diagonal
    distance[diagind(distance)].=Inf32

    # x- and y-direction of forces
    Fx = F.*(x_pairwise)./distance
    Fy = F.*(y_pairwise)./distance
    Fz = F.*(z_pairwise)./distance

    # Sum of all forces acting on each cell as a vector
    force = hcat(sum(Fx, dims = 2),sum(Fy, dims = 2),sum(Fz, dims = 2))
    
    return (force)
end 



function cell_division_3d(cell_ms, r_old, t,dt, par,simulation_type_tup)
   # xy, r, r0, r_old, r_max, t, t0

    # Choose where division will randomly occur    
    p_old = division_probability.(r_old, par["r_max"])
    p = division_probability.(cell_ms.r, par["r_max"])

    prob = (p .- p_old) ./ (1 .- p_old)
    prob[r_old .== cell_ms.r0]  = p[r_old .== cell_ms.r0]

        ####### test für dt ######
    prob = division_probability.(cell_ms.r, par["r_max"]) * dt

    random_numbers = rand(length(cell_ms.r))
    indices = findall(x -> x < 0, random_numbers - prob)


    if !isempty(indices)

        # New radius based on the area of the mother cell being two times that of the daughter cells
        r_new = cell_ms.r[indices]/2^(1/3)  #1/3 because of 3D
        # distance between the two daughter cells
        distance = rand.(Normal.((cell_ms.r[indices] - r_new)/2,0.1*(cell_ms.r[indices] - r_new)/2))

        # angle_value of cell division
        angle_value_xy = rand(length(r_new))*2*π
        # angle in z 
        angle_value_z = rand(length(r_new))*2*π

        # displacement vectors
        dx = distance .* cos.(angle_value_xy)
        dy = distance .* sin.(angle_value_xy)
        dz = distance .* sin.(angle_value_z)
    
        # displacement
        dxyz = [dx;dy;dz]
        dxyz = reshape(dxyz,(3,Int(length(dxyz)/3))) # transponiert oder nicht?? muss ich noch schauen
        xy1 = cell_ms.xy[indices,:] + dxyz'
        xy2 = cell_ms.xy[indices,:] - dxyz'
     
        # position update of mother cells
        cell_ms.xy[indices,:] = xy1
        
        # adding positions of daughter cells
        cell_ms.xy = vcat(cell_ms.xy,xy2)
        
        # radii update of mothers cells and adding radii of daughter cells
        for i in range(1,length(r_new))
            push!(cell_ms.r, r_new[i])
            push!(cell_ms.r0, r_new[i])
        end
        
        # t0 update of mothers cells and adding t0 of daughter cells
        # t0 is used to calculate cell growth
        cell_ms.t0[indices] .= t
        for i in range(1,length(cell_ms.t0[indices]),step=1)
            push!(cell_ms.t0, cell_ms.t0[indices[i]])
        end
        
        if !isempty(cell_ms.u)
            # G update of mothers cells and adding t0 of daughter cells
            # G is used to calculate cell growth
            if simulation_type_tup.tf_split
                cell_ms.v[indices] .= cell_ms.v[indices]/2
                for i in range(1,length(cell_ms.v[indices]),step=1)
                    push!(cell_ms.v, cell_ms.v[indices[i]]/2)
                end

                # N update of mothers cells and adding t0 of daughter cells
                # N is used to calculate cell growth
                cell_ms.u[indices] .=cell_ms.u[indices]/2
                for i in range(1,length(cell_ms.u[indices]),step=1)
                    push!(cell_ms.u, cell_ms.u[indices[i]]/2)
                end
            else
                for i in range(1,length(cell_ms.v[indices]),step=1)
                    push!(cell_ms.v, cell_ms.v[indices[i]])
                end

                # N update of mothers cells and adding t0 of daughter cells
                # N is used to calculate cell growth

                for i in range(1,length(cell_ms.u[indices]),step=1)
                    push!(cell_ms.u, cell_ms.u[indices[i]])
                end
            end
        end
    end
    
    return cell_ms
end


### Analysis ###

function parameter_scan(par_dict::Dict,par_name::String,par_value::Number,sim_type_tup::NamedTuple,
    cell_ms::Cell_struct,time_ss::Time_struct,data_tup::NamedTuple;path=""::String,replicates=10::Integer,ID=""::String,three_d=false::Bool)

    # create a folder for the data
    if !isdir(path*par_name)
        mkdir(path*par_name)
    end
    
    path = path*par_name*"/"

    #save sim_type_tup with save_sim_type function
    save_sim_type(sim_type_tup,"sim_type",path)

    #save par_dict with save_par function
    save_par(par_dict,"parameter",path)

    #change par_name in par_dict to par_value
    par_dict[par_name]=par_value
    #run_sim and use save_data function to save the data
    name = string(par_value)
    replace(name,"." => "_")
    for i in 1:replicates

        if sim_type_tup.randomize_initial_conditions
            cell_ms.u = par_dict["innit"] .+ rand(Normal(0,par_dict["innit_noise"]),length(cell_ms.r))
            cell_ms.v = par_dict["innit"] .+ rand(Normal(0,par_dict["innit_noise"]),length(cell_ms.r))

        end

        cell_ms_copy = deepcopy(cell_ms)
        if three_d == true
            cell_ms_copy=run_sim_3d(cell_ms_copy,time_ss,par_dict,data_tup,sim_type_tup)
        else
            cell_ms_copy=run_sim(cell_ms_copy,time_ss,par_dict,data_tup,sim_type_tup)
        end
        save_data(cell_ms_copy,par_name*"_"*name*"_rep_"*ID*string(i),path)
    end

end

function steady_state_parameter_scan(par_dict::Dict,par_name::String,par_value::Number,sim_type_tup::NamedTuple,
    cell_ms::Cell_struct,connectivity;path=""::String,dt=0.01,maxsteps=100000,
    maxerror=1e-9,replicates=10::Integer,ID=""::String,save_time_course=false::Bool)


    
    #convert the connectivity DF to a matrix
    matrix = []
    matrix = matrix_from_delaunay(connectivity, cell_ms, simulation_type_tup)

    # create a folder for the data
    if !isdir(path*par_name)
        mkdir(path*par_name)
    end
    
    path = path*par_name*"/"

    #save sim_type_tup with save_sim_type function
    save_sim_type(sim_type_tup,"sim_type",path)

    #save par_dict with save_par function
    save_par(par_dict,"parameter",path)


    #change par_name in par_dict to par_value
    par_dict[par_name]=par_value

    #run_steady_state_transcription and use save_data function to save the data
    name = string(par_value)
    replace(name,"." => "_")
    for i in 1:replicates
        
        #randomize the initial conditions

        if sim_type_tup.randomize_initial_conditions
            
            cell_ms.u = par_dict["innit"] .+ rand(Normal(0,par_dict["innit_noise"]),length(cell_ms.r))
            cell_ms.v = par_dict["innit"] .+ rand(Normal(0,par_dict["innit_noise"]),length(cell_ms.r))

        else 
            cell_ms.u = par_dict["innit"] *ones(length(cell_ms.r))
            cell_ms.v = par_dict["innit"] *ones(length(cell_ms.r))

        end

        cell_ms_copy = deepcopy(cell_ms)

        # sim to steady state and save data
        if simulation_type_tup.save_data
            cell_ms_copy, data_tup = run_steady_state(cell_ms_copy, matrix, par_dict,simulation_type_tup, dt, maxsteps, maxerror,save_time_course=true)
            save_data(cell_ms_copy,par_name*"_"*name*"_rep_"*ID*string(i),path)
            save_data(data_tup,par_name*"_"*name*"_rep_"*ID*string(i)*"_tdc",path)
        else
            cell_ms_copy = run_steady_state(cell_ms, matrix, par_dict,simulation_type_tup, dt, maxsteps, maxerror)
            save_data(cell_ms_copy,par_name*"_"*name*"_rep_"*ID*string(i),path)
        end
    end

end

function count_zeros(x)
    num_zeros = 0
    for val in x
      if val == 0.0
        num_zeros += 1
      end
    end
    return num_zeros
end

function calculate_moran_index(cell_ms, simulation_type_tup)

    # Vector with zeros
   n = length(cell_ms.u)
   x = zeros(n)

   # u is bigger than v -> 1
   x[cell_ms.u .> cell_ms.v] .= 1

   # Transform the connectivity matrix into a weight matrix
   w = matrix_from_delaunay(delaunay(cell_ms.xy'), cell_ms, simulation_type_tup)

   # y is the difference between x and the mean of x
   y = x .- mean(x)

   # numerator is the scalar product of w and y
   numerator = dot(y, w * y)

   # denominator is the sum product of y^2
   denominator = sum(y.^2)

   # Morans_I is the Moran-Index
   Morans_I= (n / sum(w)) * numerator / denominator



   num_zeros = count_zeros(x)

   #println("Number of v dominated cells: ", num_zeros)
   #println("Moran's I: ", Morans_I)
   return Morans_I
end

### Plotting ###

function plot_2d(cell_ms::Cell_struct;dims=[1,2],edges=false,radius=false,clim=(0,get_lim(cell_ms.u,cell_ms.v)),
    title="replicate: ",colorbar=true,xlabel="x",ylabel="y")


    p = scatter(cell_ms.xy[:,dims[1]],cell_ms.xy[:,dims[2]],label="",xlabel="x",ylabel="y",
    aspect_ratio=:equal,marker_z=cell_ms.u,clim = clim,colorbar = colorbar,title=title)
    if edges
        connectivity = delaunay(cell_ms.xy')
        matrix =[]
        matrix = matrix_from_delaunay(connectivity, cell_ms, simulation_type_tup)
        draw_edges_from_matrix!(matrix,cell_ms.xy)

        scatter!(cell_ms.xy[:,1],cell_ms.xy[:,2],label="",xlabel=xlabel,ylabel=ylabel,
        aspect_ratio=:equal,marker_z=cell_ms.u,clim = clim,title = title,colorbar = colorbar)
    end

    return p
end


function draw_edges_from_connectivity_3d!(connectivity,xy)

    # Packages: Combinatorics, Plots

    for k in 1:length(connectivity[1,:])    # iterate over all triangles from delaunay tesselation
        # get all edges of one triangle:
        edge_list = combinations(connectivity[:,k],2) |>collect # same as collect(combinations(connectivity[:,k],2))
        
        for i in 1:length(edge_list)    # iterate over all edges

            a = [xy[edge_list[i][1],1],xy[edge_list[i][2],1]]

            b = [xy[edge_list[i][1],2],xy[edge_list[i][2],2]]
            
            c = [xy[edge_list[i][1],3],xy[edge_list[i][2],3]]

            plot!(a,b,c,label="",alpha= 0.4,color="gray") #plot edges
        end
    end
end

function draw_edges_from_matrix_3d!(matrix,xy)

    # Packages: Plots

    # matrix = connectivity_matrix xy = cell positions

    for k in range(1,length(matrix[1,:])) # iterate over matrix
        for i in range(1,k-1)
            if matrix[k,i] == 1 # if cells i and k are connected, a line is drawn between their positions
                a = [xy[i,1],xy[k,1]]
                b = [xy[i,2],xy[k,2]]
                c = [xy[i,3],xy[k,3]]
                plot!(a,b,c,label="",alpha= 0.4,color="gray")
            end
        end
    end
end


function plot_slices(data;plot_title="plot_title",bins=6,size=(1200,800))

    # Packages: DataFrames, Plots, ColorSchemes

    test = DataFrame(data.xy, :auto)
    # create 4 equally spaced bins between the minimum and maximum values of the x2 column
    bins = range(minimum(test.x2), stop=maximum(test.x2), length=6)
    test.u = data.u

    plot_list = []   
    for i in 1:length(bins)-1
    
        filtered_test = filter(row ->  bins[i] <= row.x3 < bins[i+1], test)
        a = scatter(filtered_test.x1,filtered_test.x2,filtered_test.x3,label="",xlabel="",ylabel="",zlabel="",title="part "*string(i),aspect_ratio=:equal,marker_z=filtered_test.u,
        clim = (0, 0.1),ylim =(-4,4),xlim=(-4,4),zlim=(-4,4),colorbar = false)
        push!(plot_list,a)
    end
    push!(plot_list,scatter(test.x1,test.x2,test.x3,label="",xlabel="",ylabel="",zlabel="",title="complete",aspect_ratio=:equal,marker_z=test.u,
    clim = (0, 0.1),ylim =(-4,4),xlim=(-4,4),zlim=(-4,4),colorbar = false))
    p = plot(plot_list...,size=size,plot_title = plot_title)

    return(p)
end

function plot_time_course(x,data_matrix,time_ss;ylabel="ylabel",title_1="title_1",title_2="title_2",
    xlabel="t [au]",size_plots=(300,500),size_layout=(600,500),ylim=(0,get_lim(data_matrix)),layout_title="title layout",label="label",color="black")
    
    p = plot(ylabel=ylabel,title=title_1,ylim=ylim,size=size_plots)  
    p2 = plot(ylabel=ylabel,xlabel=xlabel,title=title_2,ylim=ylim,size=size_plots)
    counter_1 = 0
    counter_2 = 0
    for i in range(1,length(cell_ms.r),step=1)
        if x[i] == 1
            if counter_1 == 0 
                plot!(p,time_ss.t[1:length(data_matrix[:,1])],data_matrix[:,i],label=label,color=color,alpha=0.5)
                counter_1 += 1
            else
                plot!(p,time_ss.t[1:length(data_matrix[:,1])],data_matrix[:,i],label="",color=color,alpha=0.5)
            end
        else
            if counter_2 == 0 
                plot!(p2,time_ss.t[1:length(data_matrix[:,1])],data_matrix[:,i],label=label,color=color,alpha=0.5)
                counter_2 += 1
            else
                plot!(p2,time_ss.t[1:length(data_matrix[:,1])],data_matrix[:,i],label="",color=color,alpha=0.5)
            end
        end
    end
    p3 = plot(p,p2,layout=(2,1),size=size_layout,plot_title=layout_title,plot_titlefontsize=12);
    return p3
end 


function plot_time_course_diff(x,data_matrix_1,data_matrix_2,time_ss; ylabel_1="ylabel_1",ylabel_2="ylabel_2",title_1="title_1",title_2="title_2",
    xlabel_1="t [au]",xlabel_2 ="t [au]",size_plots=(300,500),size_layout=(600,500),ylim=(0,get_lim(data_matrix_1,data_matrix_2)),
    layout_title="title layout",label_1="label_1",label_2="label_2",color_1="black",color_2="grey",layout=(2,1),values=[1,0])
    
    p = plot(ylabel=ylabel_1,xlabel=xlabel_1,title=title_1,ylim=ylim,size=size_plots)  
    p2 = plot(ylabel=ylabel_2,xlabel=xlabel_2,title=title_2,ylim=ylim,size=size_plots)
    counter_1 = 0
    counter_2 = 0
    for i in range(1,length(cell_ms.r),step=1)
        if x[i] == values[1]
            if counter_1 == 0 
                plot!(p,time_ss.t[1:length(data_matrix_1[:,1])],data_matrix_1[:,i],label=label_1,color=color_1,alpha=0.5)
                counter_1 += 1
            else
                plot!(p,time_ss.t[1:length(data_matrix_1[:,1])],data_matrix_1[:,i],label="",color=color_1,alpha=0.5)
            end
        end
        if x[i] == values[2]
            if counter_2 == 0 
                plot!(p2,time_ss.t[1:length(data_matrix_2[:,1])],data_matrix_2[:,i],label=label_2,color=color_2,alpha=0.5)
                counter_2 += 1
            else
                plot!(p2,time_ss.t[1:length(data_matrix_2[:,1])],data_matrix_2[:,i],label="",color=color_2,alpha=0.5)
            end
        end
    end
    p3 = plot(p,p2,layout=layout,size=size_layout,plot_title=layout_title,plot_titlefontsize=12);
    return p3
end 



function get_lim(args...)
    max = maximum(maximum.(skipmissing.(args)))
    if  max > 1
       @warn"Maximum is > 1, returing maximum value"
       return (max) 
    elseif max > 0.1
        return (1)
    else
        return (0.1)
    end
end 


function findfirstminequal(arr,value) ### helps with some plotting of timeseries data
    for i in range(1,length(arr))
        if arr[i] >= value
            return(i)
        end
    end
end

####### time course analysis ###########

function calculate_moran_index(data, index,simulation_type_tup)

    # Vector with zeros
   n = length(data.u[index])
   x = zeros(n)

   # u is bigger than v -> 1
   x[data.u[index] .> data.v[index]] .= 1

   # Transform the connectivity matrix into a weight matrix
   w = matrix_from_delaunay(delaunay(data.xy[index]'), data, index, simulation_type_tup)

   # y is the difference between x and the mean of x
   y = x .- mean(x)

   # numerator is the scalar product of w and y
   numerator = dot(y, w * y)

   # denominator is the sum product of y^2
   denominator = sum(y.^2)

   # Morans_I is the Moran-Index
   Morans_I= (n / sum(w)) * numerator / denominator



   num_zeros = count_zeros(x)

   #println("Number of v dominated cells: ", num_zeros)
   #println("Moran's I: ", Morans_I)
   return Morans_I
end


function calculate_ratio(data,index,simulation_type_tup)

    # Vector with zeros
   n = length(data.u[index])
   x = zeros(n)

   # u is bigger than v -> 1
   x[data.u[index] .> data.v[index]] .= 1

   ratio = sum(x)/n

   return ratio
end

function calculate_ratio(data,simulation_type_tup)

    # Vector with zeros
   n = length(data.u)
   x = zeros(n)

   # u is bigger than v -> 1
   x[data.u .> data.v] .= 1

   ratio = sum(x)/n

   return ratio
end

function matrix_from_delaunay(connectivity,cell_ms,index,simulation_type_tup)
    xy = cell_ms.xy[index]'
    # Packages: Combinatorics

    m = zeros(length(xy[1,:]),length(xy[1,:])) # define graph with number of notes as cells and no edges
    for k in 1:length(connectivity[1,:])    # iterate over all triangles from delaunay tesselation
        # get all edges of one triangle:
        edge_list = combinations(connectivity[:,k],2) |>collect # same as collect(combinations(connectivity[:,k],2))

        if simulation_type_tup.check_distance # check if distance between to points is greater than the sum of their radii
            for i in 1:length(edge_list)    # iterate over all edges
                if euclidean(xy[:,edge_list[i][1]],xy[:,edge_list[i][2]]) < (cell_ms.r[index][edge_list[i][1]] + cell_ms.r[index][edge_list[i][2]])
                    m[edge_list[i][1],edge_list[i][2]] = 1    # add edges to graph
                    m[edge_list[i][2],edge_list[i][1]] = 1 
                end
            end
        else
            for i in 1:length(edge_list)    # iterate over all edges
                m[edge_list[i][1],edge_list[i][2]] = 1    # add edges to graph
                m[edge_list[i][2],edge_list[i][1]] = 1 
            end
        end
    end
return(m)
end


function get_switches(input)
    index = []
    a  = 0
    for j in range(1,length(input))
        a+=1
        if !ismissing(input[j])
            break
        end
    end      
    start = input[a]
    for i in range(a,length(input))
        if input[i] == - start
            push!(index,i)
            start = input[i]
        end
    end
    return index
end    

function calculate_switches_delta(u_array,v_array,delta)
    u_pos = false
    v_pos = false
    undecided = false
    switches = 0
    for i in range(1,length(u_array))
        if ismissing(u_array[i])
            continue
        end
        if u_array[i] > 1 - delta && v_array[i] < delta
            if v_pos == true
                switches += 1
            end
            u_pos = true
            v_pos = false

        elseif v_array[i] > 1 - delta && u_array[i] < delta

            if u_pos == true
                switches += 1
            end
            u_pos = false
            v_pos = true
        end
    end
    return switches 
end

function calculate_switches_delta_index(u_array,v_array,delta)
    u_pos = false
    v_pos = false
    undecided = false
    switches = 0
    index_array=[]
    for i in range(1,length(u_array))
        if ismissing(u_array[i])
            continue
        end
        if u_array[i] > 1 - delta && v_array[i] < delta
            if v_pos == true
                switches += 1
                push!(index_array,i)
            end
            u_pos = true
            v_pos = false

        elseif v_array[i] > 1 - delta && u_array[i] < delta

            if u_pos == true
                push!(index_array,i)
            end
            u_pos = false
            v_pos = true
        end
    end
    return switches,index_array
end


#### used for extending time series data to same length for each rep: can calculate mean and std
function fill_with_last_number(array)

    if !isnan(array[end])
        return array
    end

    for i in range(length(array),step=-1,stop=1)

        if !isnan(array[i])
            array[i:end].=array[i]
            return array
        end
    end

    end