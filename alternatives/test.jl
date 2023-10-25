using Downloads, SBMLToolkit, Graphs, SimpleWeightedGraphs, StatsPlots, InformationMeasures, GraphPlot, DataFrames, PumasQSP, ModelingToolkit, DifferentialEquations
Downloads.download("https://ftp.ebi.ac.uk/pub/databases/biomodels/repository/aaj/MODEL2107190002/5/Bakshi2020%20truncated%20minimal%20model.xml", "_model.xml")
odesys = readSBML("_model.xml", ODESystemImporter())
tspan = (0.0, 10.0)
prob = ODEProblem(odesys, [], tspan)
sol = solve(prob, saveat = 1.)
data = DataFrame(sol)
plot(sol)



# different cells with different k1 values between 0 and 1
@unpack k1 = odesys
@unpack C3b, Bb, C3, C3bB_closed, C3bBb, Factor_B, C3bB_open, Factor_D = odesys
default_value_k1 = ModelingToolkit.defaults(odesys)[k1]
n_cells = 10
k1_vals = rand(n_cells)
saveat = 0.1


function prob_func(prob, i, repeat)
    remake(prob, p = [k1 => k1_vals[i]])
end

ens_prob = EnsembleProblem(prob, prob_func = prob_func)
ens_sol = solve(ens_prob, saveat = saveat, trajectories = length(k1_vals))

plot(ens_sol)
df = DataFrame(ens_sol)


#cells/trajectories and genes/states
function diff_(a,b)
    return sum((a.-b).^2)
end
sts = [C3b, Bb, C3, C3bB_closed, C3bBb, Factor_B, C3bB_open, Factor_D]

ts = ens_sol[1].t
v = Vector{Matrix{Float64}}(undef, length(ts))
for (x, temp_t) in enumerate(ts)
    a = Matrix{Float64}(undef, length(sts), length(sts))
    for (i, s1) in enumerate(sts)    
        for (j, s2) in enumerate(sts)
            a[i,j] = diff_(ens_sol(temp_t, idxs = s1), ens_sol(select_timepoint, idxs = s2))
        end
    end
    v[x] = a
end

n_states = length(sts)
sw_graph = SimpleWeightedGraph(n_states)
nodelabel = String.(Symbol.(sts[:]))
cut_off_weight = 5
temp_t = 5
for i in 1:n_states
    for j in 1:n_states
        temp_weight  = v[temp_t][i,j]
        if temp_weight > cut_off_weight
            add_edge!(sw_graph, i, j, temp_weight)
        end
    end
end
gplot(sw_graph, nodelabel = nodelabel)