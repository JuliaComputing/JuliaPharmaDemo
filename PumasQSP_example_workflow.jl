### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 2bd1af72-2f3e-471b-bb28-32e38c1a32bd
using PumasQSP, Plots, OrdinaryDiffEq, DifferentialEquations, ModelingToolkit, DataFrames, StatsPlots, GlobalSensitivity, Statistics, OptimizationBBO;

# ╔═╡ 41f23ebc-cc05-46bc-9879-b986bb16dbc1
md"""
## *Pumas-QSP* example: The HIRES model
"""

# ╔═╡ a1e18f74-7f8f-4a05-b7d7-22b914fda4b8
md"""

This Pluto notebook is designed with the intention to demonstrate a classic PumasQSP workflow. The example model used to introduce PumasQSP is the *HIRES* model. It covers the following aspects of a QSP experiment:
1. **Model definition** -- Step 1
2. **First Model Simulations** -- Step 2
3. **Analysis in higher dimenstions**
    - Data -- Step 3
        - Generate simulated data
        - Visualise data 
    - Virtual Population analysis 
        - Generation -- Step 4
        - Test on data -- Step 5
        - Explore new conditions -- Step 6
"""

# ╔═╡ bbbfb448-771b-4b0c-b7e9-627b2fb2d2b7
md"""
## Prepare environment
"""

# ╔═╡ f960d305-eef9-46e8-93d5-b7980f70a159
md"""
For this example, we are using the following packages:
- `PumasQSP.jl`: This is used to implement the QSP workflow.
- `ModelingToolkit.jl`: This is used to describe the model.
- `OrdinaryDiffEq.jl` and `DifferentialEquations.jl`: Model simulation engine.
- `DataFrames.jl`: This is used to handle data.
- `OptimizationBBO.jl`: To define our optimization algorithm of choice.
- `Statistics.jl`: We use this to calculate the mean.
- `Plots.jl` and `StatsPlots.jl`: We use this package to visualize results.

"""

# ╔═╡ 79d9a06f-40ab-44c8-b701-7d877e70eb81
md"""
## Step 1: Model Definition
"""

# ╔═╡ 56c50a67-c7f8-4009-9f4d-e624097b1875
md"""
We define the model using the `ModelingToolkit` syntax. Here, the model is described as an `ODESystem` object which is defined by 
1. an indepenent variable time `t`,
2. the states of the systems which are `y1`, `y2`, `y3`, `y4`, `y5`, `y6`, `y7`, and `y8`,
3. the differential equations `eqs` describing the dynamics of the system.

"""

# ╔═╡ 30f6f698-7cbe-4f4c-8791-6c642c4fd816
begin 
	@variables begin
	    t 
	    y1(t) = 1.0  
	    y2(t) = 0.0
	    y3(t) = 0.0
	    y4(t) = 0.0
	    y5(t) = 0.0
	    y6(t) = 0.0
	    y7(t) = 0.0
	    y8(t) = 0.0057
	end
	@parameters begin
	    k1 = 1.71
	    k2 = 280.0
	    k3 = 8.32
	    k4 = 0.69
	    k5 = 0.43 
	    k6 = 1.81
	end
	D = Differential(t)
	eqs = [
	    D(y1) ~ (-k1*y1 + k5*y2 + k3*y3 + 0.0007),
	    D(y2) ~ (k1*y1 - 8.75*y2),
	    D(y3) ~ (-10.03*y3 + k5*y4 + 0.035*y5),
	    D(y4) ~ (k3*y2 + k1*y3 - 1.12*y4),
	    D(y5) ~ (-1.745*y5 + k5*y6 + k5*y7),
	    D(y6) ~ (-k2*y6*y8 + k4*y4 + k1*y5 - k5*y6 + k4*y7),
	    D(y7) ~ (k2*y6*y8 - k6*y7),
	    D(y8) ~ (-k2*y6*y8 + k6*y7)
	]
	@named QSPmodel = ODESystem(eqs)
end

# ╔═╡ a96e0591-0cca-4695-9a43-5f7e5204e83c
md"""
## Step 2: Simulate model
To simulate the model, we convert it to an `ODEProblem`, solve it and subsequently visualize the solution.
"""

# ╔═╡ 8ebf0e5e-9980-4511-805f-b99d14803077
begin 
	tspan = (0.0, 10.0)
	prob = ODEProblem(QSPmodel, [], tspan)
	sol = solve(prob, Rosenbrock23())
	plot(sol, grid = "off", legend = :outerright)
end

# ╔═╡ 38681cc1-e4b1-42c4-8836-79c9f145029b
md"""
## Step 3: Analysis PART 1: Data
"""

# ╔═╡ 2266bace-b5a9-40f0-85df-ed37ee4d5dd5
md"""
In this example, we consider data from the system in **three** different clinical conditions. The conditions are based on the default model setup but further specified via the given time spans and parameters.
"""

# ╔═╡ 40611a30-696e-415b-b082-1ac96eeefd0d
md"""
##### Data 1: Steady state data - Normalization/Adaption
"""

# ╔═╡ 59c3d5bf-8dc2-4f37-b7c0-1eb4d8f339ba
# Clinical condition 1
begin 
	prob1_ode = ODEProblem(QSPmodel, [], (0.0, 100000.0), [k6 => 20.0])
	prob1 = SteadyStateProblem(prob1_ode)
	sol1 = solve(prob1, DynamicSS(QNDF()))
	data1 = DataFrame(["y1(t)", "y2(t)", "y3(t)", "y4(t)", "y5(t)", "y6(t)", "y7(t)", "y8(t)"] .=> sol1)
end

# ╔═╡ 995e26f2-4aee-4239-a697-3ee895301bc8
md"""
##### Data 2: Time series data - Long term effect
"""

# ╔═╡ a3d1b87f-2246-4212-ab2d-4370e9e5278f
# Clinical condition 2
begin	
	prob2 = ODEProblem(QSPmodel, sol1.u, (0.0, 20.0), [k3 => 4.0, k4 => 0.4, k6 => 8.0])
	sol2 = solve(prob2, Rosenbrock23(), abstol=1e-9, reltol=1e-8)
	data2 = DataFrame(sol2)
end

# ╔═╡ 81beeae7-7599-45fb-877f-baf2c09c16fb
md"""
##### Data 3: Time series data -  - Short term effect
"""

# ╔═╡ 617615bd-6aa1-4994-815d-65f68d0ac474
# Clinical condition 3
begin	
	prob3 = ODEProblem(QSPmodel, sol1.u, (0.0, 10.0), [k3 => 15.0, k4 => 0.6, k6 => 2.0])
	sol3 = solve(prob3, Rosenbrock23(), abstol=1e-9, reltol=1e-8, saveat = Array(0:10))
	data3 = DataFrame(sol3)
end

# ╔═╡ 721fb15d-853b-4ee1-89ac-154b28590f53
md"""
With this data and the defined conditions, we create the `Trial` objects.
"""

# ╔═╡ 6eb29af1-29e2-4af6-baad-08d745b8f7c5
# Creating Trial for condition 1: Normalization/Adaption
begin 
	trial_1 = SteadyStateTrial(data1, QSPmodel;
                alg = DynamicSS(QNDF()),
                abstol = 1e-9, reltol = 1e-8,
                params = [k6 => 20.0],
				verbose = false
                )
end;

# ╔═╡ 0d3aae0a-0c7e-4dd9-bee0-7496a54a96d6
# Creating the Trial object for condition 2: Long term effect
begin 
	trial_2 = Trial(data2, QSPmodel;
                tspan = (0.0, 20.0),
                alg = Rosenbrock23(),
                abstol = 1e-9, reltol = 1e-8,
                forward_u0 = true,
                params = [k3 => 4.0, k4 => 0.4, k6 => 8.0],
				verbose = false
                )
end;

# ╔═╡ f6e787a8-fb1d-4802-812b-966b1a4bd334
# Creating Trial for condition 3: Short term effect
trial_3 = Trial(data3, QSPmodel;
                tspan = (0.0, 10.0),
                alg = Rosenbrock23(),
                abstol = 1e-9, reltol = 1e-8,
                forward_u0 = true,
                params = [k3 => 15.0, k4 => 0.6, k6 => 2.0],
				verbose = false
);

# ╔═╡ 6faa7f93-144b-4d1f-9aef-87679cb5025a
md"""
### Step 4. Analysis -- Generate Virtual Population 
In order to create a Virtual Population (VP), we need to define the object `InverseProblem`. For this, the following is required:
- the `QSPmodel` which is the `ODESystem` in the `ModelingToolkit` syntax,
- the `trials` which is the collection of conditions to be considered for the population `[trial_1, trial_2, trial_3`],
- the `search_space` which constrains the diversity in the populations in a semantically sensible manner. Here, this is set ot be the initial condition of the first state `y1` and the parameters `k1` and `k2`.
Once the `InverseProblem` object is created, we can implement the `vpop` call which produces a VP of the size specified in the `population_size` argument. Here, the generation is done via a stoachstic optimization procedure.
"""

# ╔═╡ c42e63de-ccec-4434-b541-23ae3b378e72
# ╠═╡ show_logs = false
begin 
	trials = SteadyStateTrials(trial_1, trial_2, trial_3);
	population_size = 10
 	inv_prob = InverseProblem(trials, QSPmodel, 
		[k1 => (0. , 5.), k2 => (200., 300.), y1 => (0. , 3.)]) 
 	vp = vpop(inv_prob, StochGlobalOpt(
		method = SingleShooting(maxiters = 750, optimizer = BBO_adaptive_de_rand_1_bin_radiuslimited())), 
		population_size = population_size)
end;

# ╔═╡ af6bf403-35a3-47dc-9312-e9c67ce6c20d
md"""
The VP generation is done once the 100% is reached. We visualise the distibution of the optimized parameters, i.e. the virtual patients below.
"""

# ╔═╡ a94c1bd6-2fbc-4412-b36d-7acc75864335
begin 
	params1 = DataFrame(vp)[!,"k1"]
	params2 = DataFrame(vp)[!,"k2"]
	params3 = DataFrame(vp)[!,"y1(t)"]
	plt1 = violin(params1, grid = :off, ylim = (0,5), label = "", ylabel = "k1", color = "purple")
	plt2 = violin(params2, grid = :off,  ylim = (200, 300),  label = "", ylabel = "k2", color = "purple")
	plt3 = violin(params3, grid = :off, label = "", ylabel = "y1", color = "purple")
	plot(plt1, plt2, plt3)
end

# ╔═╡ 09663cb0-8f4d-43cd-9d18-2a48b04f782e
md"""
### Step 5. Analysis -- Test VP on data
Now that the VP is generated, we can analyze the behavior of the system for the virtual patients -- at the population level, or individually.
First, we focus on the behavior for the conditions where we had data provided.
"""

# ╔═╡ 904c373d-a38b-4130-8b1c-dc044631b72d
md"""
There are two things to select for this:

1. Clinical condition  
   - This is to be specified by the `trial number`.
2. Format of display
   - Here, we provide different formats of display:
     - `summary mode`: this provides a summary statistic of the whole population.
     - `all trajectories`: this displays all patients as single lines.
"""

# ╔═╡ 7ca4fe49-3f21-4c30-91d4-459ff7642368
md"""
>Please selcet the trial number of the trial you want to display:
"""

# ╔═╡ 7ce39e9f-6fd7-4987-94f7-6c1da9622183
@bind trial_number html"<input type=range min=2 max=3>"

# ╔═╡ 0d95c7be-7679-4435-ab8a-0f69c59fd747
print("Trial number to be displayed: $trial_number")

# ╔═╡ 937e5c93-ffd1-4b67-8350-37d83236d082
md"""
> Please select your preferred display format:
"""

# ╔═╡ c230fc3c-17e4-4426-abfe-139311d4a26e
@bind format_type html"<input type=range min=1 max=2>"

# ╔═╡ 3faede93-38f6-432b-b39c-553bc89260ab
begin 
	formats = ["summary mode", "all trajectories"]
	temp_format = formats[format_type]
	print("Format type to be displayed: $temp_format")
	summ = true
	legend_dis = :outerright
	if format_type == 2
		summ = false
		legend_dis = :false
	end
end;

# ╔═╡ a2f66f61-bfab-4879-920c-ec120a226be7
md"""
>Do you want to visualize the data?

"""

# ╔═╡ 3f609536-e997-4dd0-9843-db7946e00feb
@bind show_data_flag html"<input type=range min=1 max=2>"

# ╔═╡ f1b7b3c6-bdae-4314-8f1e-fac63d908c6a
begin 
	show_data_flags = [false, true]
	temp_show_data_flag = show_data_flags[show_data_flag]
	print("Data shown: $temp_show_data_flag")
end;

# ╔═╡ 1d263fee-1eda-47f5-a0c0-a02e956c0a09
plot(vp, trials[trial_number]; summary = summ, show_data = temp_show_data_flag, grid = "off", ylabel = "State",legend = legend_dis)

# ╔═╡ 6752199d-2fc3-4e4e-81ba-0788ffb4a7f9
md"""
### Step 6. Exploring VP in *new conditions* 
The previous step can be seen as a form of sanity check. In the next step, we now use the VP to explore new configurations of the system, i.e. new clinical conditions. Pumas-QSP's API offers the `virtual_trial` object for this. Starting from any known clinical condition (i.e. experimentally observed condition and previously specified as a `Trial` object), we can now modify it to explore new conditions. 
"""

# ╔═╡ 37f4c431-651d-4028-95c2-3a304ecca45e
md"""
> Which trial do you want to start with?
"""

# ╔═╡ 12954f4f-ccea-4323-8b70-47ef378b4064
@bind trial_ind_new html"<input type=range min=2 max=3>"

# ╔═╡ b573d05a-b265-4f7c-b434-9a7d22c68464
print("Selected trial: $trial_ind_new")

# ╔═╡ 51194b3c-72ea-418f-9f06-66531cb50960
plot(vp, trials[trial_ind_new], legend = :outerright, grid = "off")

# ╔═╡ be61f182-d02d-41b0-8648-ee896d819ce1
md"""
Now we will create a new clinical condition to explore by creating a `virtual_trial` based on the selected trial. Available modifications are for example
- new parameters,
- new initial conditions,
- new time spans.
Combinations are also possible. For the sake of clarity, we are sticking with individual changes in this workflow.
"""

# ╔═╡ 70819d97-7187-4b92-87b7-29ae202fa0f6
md"""
#### 6.1. Explore new parameters
For example, here we change the value of parameter k1 where the search space of previous virtual population was set to the range zero and five.
"""


# ╔═╡ 07d35ca9-40b3-43f5-8bd1-eb4c06710a3e
md"""
> What should the new value for parameter k1 be? Selcect between zero and five.
"""

# ╔═╡ 79ecafd1-745c-4889-a1b7-d9eaceb9a5d7
@bind param_val_new html"<input type=range min=0 max=5>"

# ╔═╡ 88592d7f-8d41-4943-b381-36f9aedb9b8f
print("Selected parameter set to: $param_val_new")

# ╔═╡ 7f5842a5-c184-4f59-92cd-fb616e546462
md"""
In the following we visualize the behaviour for the whole population and in summary mode.
"""

# ╔═╡ 7d9469a1-cd1d-4696-9661-f3ad53cb72fa
begin
	new_param = [ k1 => param_val_new]
	vt_1 = virtual_trial(trials[trial_ind_new], QSPmodel; params = new_param)
	plot(vp, vt_1, grid = "off", legend = :outerright)
end

# ╔═╡ 4d8631a5-b770-4659-b7f5-7d37511d730f

md"""
#### 6.2. Explore new initial condition
Here you can select, which initial condition to vary and then to which new value it should be set.
"""

# ╔═╡ 1362250a-c20a-4709-9c2c-55bdbf8b3e9d
md"""
> Which initial condition do you want to vary?
"""

# ╔═╡ 398ad035-21c7-4055-a206-806c69bcaa1a
@bind ic_ind_new html"<input type=range min=2 max=8>"

# ╔═╡ d91ec775-e2ef-4b78-96e3-eeb789ca7f55
print("Selected initial condition to be varied: $ic_ind_new")

# ╔═╡ d4837b7c-2d4f-4289-9f06-43922856fb1d
md"""
> What should the new value for the selected initial condition be?
"""

# ╔═╡ 34d5c411-ce58-4d7c-9632-7f4b91a1d480
@bind ic_val_new html"<input type=range min=0 max=5>"

# ╔═╡ a0b29e6b-5f52-4e75-85bb-cf3a6e0b2d20
print("Selected initial condition set to: $ic_val_new")

# ╔═╡ 60a1677c-ffb2-4b73-97a0-a3cd2c2a5580
begin
	new_ic = [ y2 => ic_val_new]
	if ic_ind_new == 3
		new_ic = [ y3 => ic_val_new]
	elseif ic_ind_new == 4
		new_ic = [ y4 => ic_val_new]
	elseif ic_ind_new == 5
		new_ic = [ y5 => ic_val_new]
	elseif ic_ind_new == 6
		new_ic = [ y6 => ic_val_new]
	elseif ic_ind_new == 7
		new_ic = [ y7 => ic_val_new]
	elseif ic_ind_new == 8
		new_ic = [ y8 => ic_val_new]
	end
	vt_2 = virtual_trial(trials[trial_ind_new], QSPmodel; u0 = new_ic)
	plot(vp, vt_2, grid = "off", legend = :outerright)
end

# ╔═╡ f580f2f6-60a0-46e3-b8ff-57090a607a03

md"""
#### 6.3. Explore new time span
Here you can set a new time span.
"""

# ╔═╡ 6ba09f38-313d-4e90-a83b-f902f42a4ea4
md"""
> Please specify you end point of the time span:
"""

# ╔═╡ 265dff45-b8e6-4e28-aac7-b557ff52eddb
@bind t_end html"<input type=range min=1 max=100>"

# ╔═╡ 75fbd911-9995-48ce-b61c-2f9083b5aa16
print("Selected end point of simulation: $t_end")

# ╔═╡ 0b376a7a-7906-438c-9738-56390a03af54
begin
	new_tspan = (0, t_end)
	vt_3 = virtual_trial(trials[trial_ind_new], QSPmodel; 
	tspan = new_tspan, saveat = 1:t_end)
	plot(vp, vt_3, grid = "off", legend = :outerright)
end

# ╔═╡ eeb8a453-e9af-400e-a7c1-5c994a97b2f8
md"""
#### 6.4 Manual trial definition
Here is one more example, were you can see the `virtual_trial` object and change several things at the same time.
"""

# ╔═╡ eea3747a-23b9-4653-a789-2ab7f6f3055f
# Creating virtual trial
vt_4 = virtual_trial(trials[trial_ind_new], QSPmodel; 
	u0 =[y3 => 2.1], 
	params = [k1 => 1.7], 
	tspan = (0, 100.),
	saveat = 1:100);

# ╔═╡ b74abd28-7fd0-4e01-b821-c53acef1f715
plot(vp, vt_4, grid = "off", legend = :outerright)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
GlobalSensitivity = "af5da776-676b-467e-8baf-acd8249e4f0f"
ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
OptimizationBBO = "3e6eede4-6085-4f62-9a71-46d9bc1eb92b"
OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PumasQSP = "5526fe93-3145-4cc2-889a-1fb5ef00f63e"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"

[compat]
DataFrames = "~1.6.1"
DifferentialEquations = "~7.10.0"
GlobalSensitivity = "~2.2.1"
ModelingToolkit = "~8.70.0"
OptimizationBBO = "~0.1.5"
OrdinaryDiffEq = "~6.57.0"
Plots = "~1.39.0"
PumasQSP = "~2.4.5"
StatsPlots = "~0.15.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "e69a8cfa00385f2d22022cd613b03d8d667be783"

[[deps.ADTypes]]
git-tree-sha1 = "5d2e21d7b0d8c22f67483ef95ebdc39c0e6b6003"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.4"

[[deps.AMD]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "d4b99dd70d7136fe75ec74d072191d688448d39c"
uuid = "14f7f29c-3bd6-536c-9a0b-7339e30b5a3e"
version = "0.5.2"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "1db808640788803779048afdbe359a172e62927e"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.32.4"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractMCMC]]
deps = ["BangBang", "ConsoleProgressMonitor", "Distributed", "LogDensityProblems", "Logging", "LoggingExtras", "ProgressLogging", "Random", "StatsBase", "TerminalLoggers", "Transducers"]
git-tree-sha1 = "87e63dcb990029346b091b170252f3c416568afc"
uuid = "80f14c24-f653-4e6a-9b94-39d6b0f70001"
version = "4.4.2"

[[deps.AbstractPPL]]
deps = ["AbstractMCMC", "DensityInterface", "Random", "Setfield", "SparseArrays"]
git-tree-sha1 = "caa9b62583577b0d6b222f11f54aa29fabbdb5ca"
uuid = "7a57a42e-76ec-4ea3-a279-07e840d6d9cf"
version = "0.6.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdvancedHMC]]
deps = ["AbstractMCMC", "ArgCheck", "DocStringExtensions", "InplaceOps", "LinearAlgebra", "LogDensityProblems", "LogDensityProblemsAD", "ProgressMeter", "Random", "Requires", "Setfield", "SimpleUnPack", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "acbe805c3078ba0057bb56985248bd66bce016b1"
uuid = "0bf59076-c3b1-5ca4-86bd-e02cd72cde3d"
version = "0.5.5"

    [deps.AdvancedHMC.extensions]
    AdvancedHMCCUDAExt = "CUDA"
    AdvancedHMCMCMCChainsExt = "MCMCChains"
    AdvancedHMCOrdinaryDiffEqExt = "OrdinaryDiffEq"

    [deps.AdvancedHMC.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
    OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"

[[deps.AdvancedMH]]
deps = ["AbstractMCMC", "Distributions", "FillArrays", "LinearAlgebra", "LogDensityProblems", "Random", "Requires"]
git-tree-sha1 = "b2a1602952739e589cf5e2daff1274a49f22c9a4"
uuid = "5b7e9947-ddc0-4b3f-9b55-0d8042f74170"
version = "0.7.5"
weakdeps = ["DiffResults", "ForwardDiff", "MCMCChains", "StructArrays"]

    [deps.AdvancedMH.extensions]
    AdvancedMHForwardDiffExt = ["DiffResults", "ForwardDiff"]
    AdvancedMHMCMCChainsExt = "MCMCChains"
    AdvancedMHStructArraysExt = "StructArrays"

[[deps.AdvancedPS]]
deps = ["AbstractMCMC", "Distributions", "Libtask", "Random", "Random123", "StatsFuns"]
git-tree-sha1 = "4d73400b3583147b1b639794696c78202a226584"
uuid = "576499cb-2369-40b2-a588-c64705576edc"
version = "0.4.3"

[[deps.AdvancedVI]]
deps = ["Bijectors", "Distributions", "DistributionsAD", "DocStringExtensions", "ForwardDiff", "LinearAlgebra", "ProgressMeter", "Random", "Requires", "StatsBase", "StatsFuns", "Tracker"]
git-tree-sha1 = "1f919a9c59cf3dfc68b64c22c453a2e356fca473"
uuid = "b5ca4192-6429-45e5-a2d9-87aec30a685c"
version = "0.2.4"

[[deps.AlgebraicMultigrid]]
deps = ["CommonSolve", "LinearAlgebra", "Printf", "Reexport", "SparseArrays"]
git-tree-sha1 = "796eedcb42226861a51d92d28ee82d4985ee860b"
uuid = "2169fc97-5a83-5252-b627-83903c6c433c"
version = "0.5.1"

[[deps.AlphaStableDistributions]]
deps = ["Distributions", "FileIO", "Interpolations", "JLD2", "LinearAlgebra", "Random", "SpecialFunctions", "StaticArrays", "Statistics", "StatsBase", "ToeplitzMatrices"]
git-tree-sha1 = "ad5789d4d1a7e38a9fb1e5e02369eeb799833b2e"
uuid = "f20549b4-2d50-407f-863c-cdd202ba59a3"
version = "1.1.6"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f83ec24f76d4c8f525099b2ac475fc098138ec31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.11"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "0d61921af2799487b80453a44abb57db7a0c1381"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.4.1"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.ArraysOfArrays]]
deps = ["Adapt", "ChainRulesCore", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "c59b725b0aadf7df93fb3de05b5e1b14029af2da"
uuid = "65a8f2f4-9b39-5baf-92e2-a9cc46fdf018"
version = "0.6.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "c06a868224ecba914baa6942988e2f2aade419be"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "0.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "0b816941273b5b162be122a6c94d706e3b3125ca"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.38"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables"]
git-tree-sha1 = "e28912ce94077686443433c2800104b061a827ed"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.39"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijections]]
git-tree-sha1 = "71281c0c28f97e0adeed24fdaa6bf7d37177f297"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.5"

[[deps.Bijectors]]
deps = ["ArgCheck", "ChainRules", "ChainRulesCore", "ChangesOfVariables", "Compat", "Distributions", "Functors", "InverseFunctions", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "MappedArrays", "Random", "Reexport", "Requires", "Roots", "SparseArrays", "Statistics"]
git-tree-sha1 = "8eacff457e5b8c13a97848484ad650dabbffa0fc"
uuid = "76274a88-744f-5084-9051-94815aaf08c4"
version = "0.13.7"

    [deps.Bijectors.extensions]
    BijectorsDistributionsADExt = "DistributionsAD"
    BijectorsForwardDiffExt = "ForwardDiff"
    BijectorsLazyArraysExt = "LazyArrays"
    BijectorsReverseDiffExt = "ReverseDiff"
    BijectorsTrackerExt = "Tracker"
    BijectorsZygoteExt = "Zygote"

    [deps.Bijectors.weakdeps]
    DistributionsAD = "ced4e74d-a319-5a8a-b0ac-84af2272839c"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.BlackBoxOptim]]
deps = ["CPUTime", "Compat", "Distributed", "Distributions", "HTTP", "JSON", "LinearAlgebra", "Printf", "Random", "SpatialIndexing", "StatsBase"]
git-tree-sha1 = "136079f37e3514ec691926093924b591a8842f5d"
uuid = "a134a8b2-14d6-55f6-9291-3336d3ab0209"
version = "0.6.2"

[[deps.BoundaryValueDiffEq]]
deps = ["ArrayInterface", "BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NonlinearSolve", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "TruncatedStacktraces", "UnPack"]
git-tree-sha1 = "f7392ce20e6dafa8fee406142b1764de7d7cd911"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "4.0.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "601f7e7b3d36f18790e2caf83a882d88e9b71ff1"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.4"

[[deps.CPUTime]]
git-tree-sha1 = "2dcc50ea6a0a1ef6440d6eecd0fe3813e5671f45"
uuid = "a9c8d775-2e2e-55fc-8582-045d282d599e"
version = "1.0.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "32abd86e3c2025db5172aa182b982debed519834"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.1"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSTParser]]
deps = ["Tokenize"]
git-tree-sha1 = "3ddd48d200eb8ddf9cb3e0189fc059fd49b97c1f"
uuid = "00ebfdb7-1f24-5e51-bd34-a7502290713f"
version = "3.3.6"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "44dbf560808d49041989b8a96cae4cffbeb7966a"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.11"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.Catalyst]]
deps = ["DataStructures", "DiffEqBase", "DocStringExtensions", "Graphs", "JumpProcesses", "LaTeXStrings", "Latexify", "LinearAlgebra", "MacroTools", "ModelingToolkit", "Parameters", "Reexport", "Requires", "SparseArrays", "SymbolicUtils", "Symbolics", "Unitful"]
git-tree-sha1 = "e38c5c3c9522f31e3089ffb4f022edda251641c7"
uuid = "479239e8-5488-4da2-87a7-35f2df7eef83"
version = "13.4.1"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "SparseInverseSubset", "Statistics", "StructArrays", "SuiteSparse"]
git-tree-sha1 = "dbeca245b0680f5393b4e6c40dcead7230ab0b3b"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.54.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "2fba81a302a7be671aefe194f0525ef231104e7f"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.8"
weakdeps = ["InverseFunctions"]

    [deps.ChangesOfVariables.extensions]
    ChangesOfVariablesInverseFunctionsExt = "InverseFunctions"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "b86ac2c5543660d238957dbde5ac04520ae977a7"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "02aa26a4cf76381be7f66e020a3eddeb27b0a092"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "PrecompileTools", "URIs"]
git-tree-sha1 = "532c4185d3c9037c0237546d817858b23cf9e071"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.12"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ComponentArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "Functors", "LinearAlgebra", "PackageExtensionCompat", "StaticArrayInterface", "StaticArraysCore"]
git-tree-sha1 = "3ee12062466a73b5f9c04113375e4934f174bc2c"
uuid = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
version = "0.15.3"

    [deps.ComponentArrays.extensions]
    ComponentArraysAdaptExt = "Adapt"
    ComponentArraysConstructionBaseExt = "ConstructionBase"
    ComponentArraysGPUArraysExt = "GPUArrays"
    ComponentArraysRecursiveArrayToolsExt = "RecursiveArrayTools"
    ComponentArraysReverseDiffExt = "ReverseDiff"
    ComponentArraysSciMLBaseExt = "SciMLBase"
    ComponentArraysTrackerExt = "Tracker"
    ComponentArraysZygoteExt = "Zygote"

    [deps.ComponentArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    GPUArrays = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SciMLBase = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "5372dbbf8f0bdb8c700db5367132925c0771ef7e"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.1"

[[deps.ConsoleProgressMonitor]]
deps = ["Logging", "ProgressMeter"]
git-tree-sha1 = "3ab7b2136722890b9af903859afcf457fa3059e8"
uuid = "88cd18e8-d9cc-4ea6-8889-5259c0d15c8b"
version = "0.1.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.ContextVariablesX]]
deps = ["Compat", "Logging", "UUIDs"]
git-tree-sha1 = "25cc3803f1030ab855e383129dcd3dc294e322cc"
uuid = "6add18c4-b38d-439d-96f6-d6bc489c04c5"
version = "0.1.3"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Copulas]]
deps = ["AlphaStableDistributions", "Cubature", "Distributions", "ForwardDiff", "GSL", "MvNormalCDF", "Random", "Roots", "SpecialFunctions", "StatsBase", "TaylorSeries"]
git-tree-sha1 = "bf5b74cd9b240eea4b713bf3b04211dabf827688"
uuid = "ae264745-0b69-425e-9d9d-cf662c5eec93"
version = "0.1.10"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.Cubature]]
deps = ["Cubature_jll"]
git-tree-sha1 = "c3f4b3b38abd7b5c3ccf59adab2568212e7530d3"
uuid = "667455a9-e2ce-5579-9412-b964f529a492"
version = "1.5.1"

[[deps.Cubature_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0fe9efb84e3eb7b14f885a95aaa0ed50c7e839c8"
uuid = "7bc98958-0e37-5d67-a6ac-a3a19030071a"
version = "1.0.5+0"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "f7f4319567fe769debfcf7f8c03d8da1dd4e2fb0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.9"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataDrivenDiffEq]]
deps = ["CommonSolve", "DataInterpolations", "DiffEqBase", "DocStringExtensions", "LinearAlgebra", "MLUtils", "ModelingToolkit", "Parameters", "ProgressMeter", "QuadGK", "Random", "RecipesBase", "Reexport", "Setfield", "Statistics", "StatsBase", "SymbolicUtils", "Symbolics"]
git-tree-sha1 = "912924cc547e80f65e1a581d68e4a604e6b32bd3"
uuid = "2445eb08-9709-466a-b3fc-47e12bd697a2"
version = "1.3.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataInterpolations]]
deps = ["LinearAlgebra", "RecipesBase", "RecursiveArrayTools", "Reexport", "Requires"]
git-tree-sha1 = "e868551161bcf933c97b3a985e25ce32e7dd04e5"
uuid = "82cc6244-b520-54b8-b5a6-8a565e85f1d0"
version = "4.2.0"

    [deps.DataInterpolations.extensions]
    DataInterpolationsChainRulesCoreExt = "ChainRulesCore"
    DataInterpolationsOptimExt = "Optim"
    DataInterpolationsRegularizationToolsExt = "RegularizationTools"
    DataInterpolationsSymbolicsExt = "Symbolics"

    [deps.DataInterpolations.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Optim = "429524aa-4258-5aef-a3af-852621145aeb"
    RegularizationTools = "29dad682-9a27-4bc3-9c72-016788665182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[[deps.DataSets]]
deps = ["AbstractTrees", "Base64", "Markdown", "REPL", "ReplMaker", "ResourceContexts", "SHA", "TOML", "UUIDs"]
git-tree-sha1 = "bbb3ea880461e62709e8cb181856de19c383e99c"
uuid = "c9661210-8a83-48f0-b833-72e62abce419"
version = "0.2.9"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "92f9967824241fe37d02e6dd76afb73b562a3fcf"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.43.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DataStructures", "DocStringExtensions", "EnumX", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "Requires", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "0d9982e8dee851d519145857e79a17ee33ede154"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.130.0"

    [deps.DiffEqBase.extensions]
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"
    DiffEqBaseZygoteExt = "Zygote"

    [deps.DiffEqBase.weakdeps]
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "Functors", "LinearAlgebra", "Markdown", "NLsolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "acc53f895588767cbb296d3d8581ebd203524a2e"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.33.0"
weakdeps = ["OrdinaryDiffEq", "Sundials"]

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "57ed4597a309c5b2a10cab5f9813adcb78f92117"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.19.0"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiableFlatten]]
deps = ["ChainRulesCore", "LinearAlgebra", "NamedTupleTools", "OrderedCollections", "Requires", "SparseArrays"]
git-tree-sha1 = "f4dc2c1d994c7e2e602692a7dadd2ac79212c3a9"
uuid = "c78775a3-ee38-4681-b694-0504db4f5dc7"
version = "0.1.1"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "96a19f498504e4a3b39524196b73eb60ccef30e9"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.10.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "b6def76ffad15143924a2199f72a5cd883a2e8a9"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.9"
weakdeps = ["SparseArrays"]

    [deps.Distances.extensions]
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "3d5873f811f582873bb9871fc9c451784d5dc8c7"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.102"
weakdeps = ["ChainRulesCore", "DensityInterface"]

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

[[deps.DistributionsAD]]
deps = ["Adapt", "ChainRules", "ChainRulesCore", "Compat", "Distributions", "FillArrays", "LinearAlgebra", "PDMats", "Random", "Requires", "SpecialFunctions", "StaticArrays", "StatsFuns", "ZygoteRules"]
git-tree-sha1 = "d61f08c7bd15c5ab215fd7a2eb61c1ae15d8ff5e"
uuid = "ced4e74d-a319-5a8a-b0ac-84af2272839c"
version = "0.6.53"

    [deps.DistributionsAD.extensions]
    DistributionsADForwardDiffExt = "ForwardDiff"
    DistributionsADLazyArraysExt = "LazyArrays"
    DistributionsADReverseDiffExt = "ReverseDiff"
    DistributionsADTrackerExt = "Tracker"

    [deps.DistributionsAD.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "51b4b84d33ec5e0955b55ff4b748b99ce2c3faa9"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.7"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPPL]]
deps = ["AbstractMCMC", "AbstractPPL", "BangBang", "Bijectors", "ChainRulesCore", "ConstructionBase", "Distributions", "DocStringExtensions", "LinearAlgebra", "LogDensityProblems", "MacroTools", "OrderedCollections", "Random", "Requires", "Setfield", "Test", "ZygoteRules"]
git-tree-sha1 = "e18dc6b49a5cfabe26aa018cf44c52b7ecd16cf1"
uuid = "366bfd00-2699-11ea-058f-f148b4cae6d8"
version = "0.23.18"
weakdeps = ["MCMCChains"]

    [deps.DynamicPPL.extensions]
    DynamicPPLMCMCChainsExt = ["MCMCChains"]

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "fea68c84ba262b121754539e6ea0546146515d4f"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.3"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EllipticalSliceSampling]]
deps = ["AbstractMCMC", "ArrayInterface", "Distributions", "Random", "Statistics"]
git-tree-sha1 = "973b4927d112559dc737f55d6bf06503a5b3fc14"
uuid = "cad2338a-1db2-11e9-3401-43bc07c9ede2"
version = "1.1.0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EnzymeCore]]
deps = ["Adapt"]
git-tree-sha1 = "3250001c57b9a3e18e1e5a257fb9ec2c012286c6"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.5.3"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ErrorfreeArithmetic]]
git-tree-sha1 = "d6863c556f1142a061532e79f611aa46be201686"
uuid = "90fa49ef-747e-5e6f-a989-263ba693cf1a"
version = "0.5.2"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "602e4585bcbd5a25bc06f514724593d13ff9e862"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.25.0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "0fa3b52a04a4e210aeb1626def9c90df3ae65268"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.1.0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "b4fbdd20c889804969571cc589900803edda16b7"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FLoops]]
deps = ["BangBang", "Compat", "FLoopsBase", "InitialValues", "JuliaVariables", "MLStyle", "Serialization", "Setfield", "Transducers"]
git-tree-sha1 = "ffb97765602e3cbe59a0589d237bf07f245a8576"
uuid = "cc61a311-1640-44b5-9fba-1b764f453329"
version = "0.2.1"

[[deps.FLoopsBase]]
deps = ["ContextVariablesX"]
git-tree-sha1 = "656f7a6859be8673bf1f35da5670246b923964f7"
uuid = "b9860ae5-e623-471e-878b-f6a53c775ea6"
version = "0.1.1"

[[deps.FMI]]
deps = ["ChainRulesCore", "DiffEqCallbacks", "DifferentialEquations", "Downloads", "FMIExport", "FMIImport", "ForwardDiff", "ProgressMeter", "Requires", "SciMLBase"]
git-tree-sha1 = "0f1c10044f7ded752b6b91b4e390e3365bb88b9c"
uuid = "14a09403-18e3-468f-ad8a-74f8dda2d9ac"
version = "0.10.2"

[[deps.FMICore]]
git-tree-sha1 = "9ae3da46a35de556cb5c9b574aa801d1c7b0ca05"
uuid = "8af89139-c281-408e-bce2-3005eb87462f"
version = "0.12.1"

[[deps.FMIExport]]
deps = ["Dates", "EzXML", "FMICore", "UUIDs"]
git-tree-sha1 = "1c786f943a3bd3d355f36bf2b8b50617d78b9692"
uuid = "31b88311-cab6-44ed-ba9c-fe5a9abbd67a"
version = "0.1.9"

[[deps.FMIImport]]
deps = ["ChainRulesCore", "Downloads", "EzXML", "FMICore", "ForwardDiff", "ForwardDiffChainRules", "Libdl", "ZipFile"]
git-tree-sha1 = "bff49e94cf2a822e92d0cbb032863589bdb47681"
uuid = "9fcbc62e-52a0-44e9-a616-1359a0008194"
version = "0.12.2"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "aa9925a229d45fe3018715238956766fa21804d1"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.6"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "b12f05108e405dadcc2aff0008db7f831374e051"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.0"

[[deps.FastRounding]]
deps = ["ErrorfreeArithmetic", "LinearAlgebra"]
git-tree-sha1 = "6344aa18f654196be82e62816935225b3b9abe44"
uuid = "fa42c844-2597-5d31-933b-ebd51ab2693f"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "9f00e42f8d99fdde64d40c8ea5d14269a2e2c1aa"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.21"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "a20eaa3ad64254c61eeb5f230d9306e937405434"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.6.1"
weakdeps = ["SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "c6e4a1fbe73b31a3dea94b1da449503b8830c306"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.21.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.ForwardDiffChainRules]]
deps = ["ChainRulesCore", "DifferentiableFlatten", "ForwardDiff", "MacroTools"]
git-tree-sha1 = "45de5ee8818983436ca91f9f0bebfe405d4fdfef"
uuid = "c9556dd2-1aed-4cfe-8560-1557cf593001"
version = "0.1.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9a68d75d466ccc1218d0552a8e1631151c569545"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.5"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.GSL]]
deps = ["GSL_jll", "Libdl", "Markdown"]
git-tree-sha1 = "3ebd07d519f5ec318d5bc1b4971e2472e14bd1f0"
uuid = "92c85e6c-cbff-5e0c-80f7-495c94daaecd"
version = "1.0.1"

[[deps.GSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "56f1e2c9e083e0bb7cf9a7055c280beb08a924c0"
uuid = "1b77fbbe-d8ee-58f0-85f9-836ddc23a7a4"
version = "2.7.2+0"

[[deps.GenericLinearAlgebra]]
deps = ["LinearAlgebra", "Printf", "Random", "libblastrampoline_jll"]
git-tree-sha1 = "02be7066f936af6b04669f7c370a31af9036c440"
uuid = "14197337-ba66-59df-a3e3-ca00e7dcff7a"
version = "0.3.11"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "d53480c0793b13341c40199190f92c611aa2e93c"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.2"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "424a5a6ce7c5d97cca7bcc4eac551b97294c54af"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.9"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Glob]]
git-tree-sha1 = "97285bbd5230dd766e9ef6749b80fc617126d496"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.1"

[[deps.GlobalSensitivity]]
deps = ["Combinatorics", "Copulas", "Distributions", "FFTW", "ForwardDiff", "KernelDensity", "LinearAlgebra", "Parameters", "QuasiMonteCarlo", "Random", "RecursiveArrayTools", "Statistics", "StatsBase", "ThreadsX", "Trapz"]
git-tree-sha1 = "e5484bc811d0335814e3f7d3815dcea19d261fbc"
uuid = "af5da776-676b-467e-8baf-acd8249e4f0f"
version = "2.2.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "899050ace26649433ef1af25bc17a815b3db52b7"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.9.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "ExprTools", "Logging", "MultivariatePolynomials", "Primes", "Random", "SIMD", "SnoopPrecompile"]
git-tree-sha1 = "44f595de4f6485ab5ba71fe257b5eadaa3cf161e"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.4.4"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5eab648309e2e060198b45820af1a37182de3cce"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "eb8fed28f4994600e29beef49744639d985a04b2"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.16"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.HypothesisTests]]
deps = ["Combinatorics", "Distributions", "LinearAlgebra", "Printf", "Random", "Rmath", "Roots", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "4b5d5ba51f5f473737ed9de6d8a7aa190ad8c72f"
uuid = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
version = "0.11.0"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InplaceOps]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "50b41d59e7164ab6fda65e71049fee9d890731ff"
uuid = "505f98c9-085e-5b2c-8e89-488be7bf1f34"
version = "0.3.0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ad37c091f7d7daf900963171600d7c1c5c3ede32"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "FastRounding", "LinearAlgebra", "Markdown", "Random", "RecipesBase", "RoundingEmulator", "SetRounding", "StaticArrays"]
git-tree-sha1 = "5ab7744289be503d76a944784bac3f2df7b809af"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.20.9"

[[deps.IntervalSets]]
deps = ["Dates", "Random"]
git-tree-sha1 = "8e59ea773deee525c99a8018409f64f19fb719e6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.7"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "68772f49f54b479fa88ace904f6127f0a3bb2e46"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.12"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "4ced6667f9974fc5c5943fa5e2ef1ca43ea9e450"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.8.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "Requires", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "c11d691a0dc8e90acfa4740d293ade57f68bfdbb"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.35"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.JuliaFormatter]]
deps = ["CSTParser", "CommonMark", "DataStructures", "Glob", "Pkg", "PrecompileTools", "Tokenize"]
git-tree-sha1 = "c0ac2f93058fd148267154c42dd7dc79869f7923"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "1.0.37"

[[deps.JuliaSimModelOptimizer]]
deps = ["ArraysOfArrays", "CSV", "CommonSolve", "ComponentArrays", "DataDrivenDiffEq", "DataFrames", "DataInterpolations", "DataSets", "Dates", "DiffEqBase", "DiffEqCallbacks", "DifferentialEquations", "Distances", "Distributed", "Distributions", "DistributionsAD", "FMI", "FMIImport", "ForwardDiff", "GenericLinearAlgebra", "GeometryBasics", "GlobalSensitivity", "HypothesisTests", "LabelledArrays", "LazyArrays", "LineSearches", "LinearAlgebra", "Logging", "Lux", "MCMCChains", "ModelingToolkit", "MultiKDE", "NNlib", "NoiseRobustDifferentiation", "NonlinearSolve", "Optimization", "OptimizationBBO", "OptimizationOptimJL", "OrdinaryDiffEq", "Preferences", "PrettyTables", "Printf", "ProgressLogging", "QuasiMonteCarlo", "Random", "RecipesBase", "RecursiveArrayTools", "RuntimeGeneratedFunctions", "SBML", "SBMLToolkit", "SciMLBase", "SciMLNLSolve", "Setfield", "SpecialFunctions", "Statistics", "StatsBase", "SteadyStateDiffEq", "StructArrays", "SymbolicUtils", "Symbolics", "Tables", "Turing", "UUIDs", "YAML"]
git-tree-sha1 = "fe518646862dc8eba47bb35feeb5c7ddf0ed02df"
uuid = "e76a656c-db09-4a9a-8f5e-6a4d9ca4a950"
version = "6.2.0"

[[deps.JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "3de1d557e382cad270d921fbc22351f5628e7b1f"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.8.0"
weakdeps = ["FastBroadcast"]

    [deps.JumpProcesses.extensions]
    JumpProcessFastBroadcastExt = "FastBroadcast"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "884c2968c2e8e7e6bf5956af88cb46aa745c854b"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.1"

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "LinearAlgebra", "MacroTools", "PrecompileTools", "Requires", "SparseArrays", "StaticArrays", "UUIDs", "UnsafeAtomics", "UnsafeAtomicsLLVM"]
git-tree-sha1 = "4c5875e4c228247e1c2b087669846941fb6e0118"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.8"
weakdeps = ["EnzymeCore"]

    [deps.KernelAbstractions.extensions]
    EnzymeExt = "EnzymeCore"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "90442c50e202a5cdf21a7899c66b240fdef14035"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.7"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "17e462054b42dcdda73e9a9ba0c67754170c88ae"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.4"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "4ea2928a96acfcf8589e6cd1429eff2a3a82c366"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "6.3.0"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "e7c01b69bcbcb93fd4cbc3d0fea7d229541e18d2"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.26+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LRUCache]]
git-tree-sha1 = "d36130483e3b6e4cd88d81633b596563264f15db"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.5.0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "cd04158424635efd05ff38d5f55843397b7416a9"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.14.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LatinHypercubeSampling]]
deps = ["Random", "StableRNGs", "StatsBase", "Test"]
git-tree-sha1 = "825289d43c753c7f1bf9bed334c253e9913997f8"
uuid = "a5e1c1ea-c99a-51d3-a14d-a9a37257b02d"
version = "1.9.0"

[[deps.LatticeRules]]
deps = ["Random"]
git-tree-sha1 = "7f5b02258a3ca0221a6a9710b0a0a2e8fb4957fe"
uuid = "73f95e8e-ec14-4e6a-8b18-0d2e271c4e55"
version = "0.0.1"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "88b8f66b604da079a627b6fb2860d3704a6729a1"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.14"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "MatrixFactorizations", "SparseArrays"]
git-tree-sha1 = "80e1594bc9e58e0c889c0cc3b9e4374cee9ec9b1"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "1.8.0"
weakdeps = ["StaticArrays"]

    [deps.LazyArrays.extensions]
    LazyArraysStaticArraysExt = "StaticArrays"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "fb6803dafae4a5d62ea5cab204b1e657d9737e7f"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.2.0"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtask]]
deps = ["FunctionWrappers", "LRUCache", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "345a40c746404dd9cb1bbc368715856838ab96f2"
uuid = "6f1fad26-d15e-5dc8-ae53-837a1d7b8c9f"
version = "0.8.6"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LimitedLDLFactorizations]]
deps = ["AMD", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "422ddab7cac3f51acc5ae23c96e710d07df9a73e"
uuid = "f5a24dde-3ab7-510b-b81b-6a72c6098d3b"
version = "0.5.0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearMaps]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6698ab5e662b47ffc63a82b2f43c1cee015cf80d"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.11.0"
weakdeps = ["ChainRulesCore", "SparseArrays", "Statistics"]

    [deps.LinearMaps.extensions]
    LinearMapsChainRulesCoreExt = "ChainRulesCore"
    LinearMapsSparseArraysExt = "SparseArrays"
    LinearMapsStatisticsExt = "Statistics"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ConcreteStructs", "DocStringExtensions", "EnumX", "EnzymeCore", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "Libdl", "LinearAlgebra", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "ba01f7a97d3d8bd711b2c00a8a68c887d8a85c9d"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.8.1"

    [deps.LinearSolve.extensions]
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveEnzymeExt = "Enzyme"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMKLExt = "MKL_jll"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = "Pardiso"

    [deps.LinearSolve.weakdeps]
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    MKL_jll = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"

[[deps.LogDensityProblems]]
deps = ["ArgCheck", "DocStringExtensions", "Random"]
git-tree-sha1 = "f9a11237204bc137617194d79d813069838fcf61"
uuid = "6fdf6af0-433a-55f7-b3ed-c6c6e0b8df7c"
version = "2.1.1"

[[deps.LogDensityProblemsAD]]
deps = ["DocStringExtensions", "LogDensityProblems", "Requires", "SimpleUnPack"]
git-tree-sha1 = "a0512ad65f849536b5a52e59b05c59c25cdad943"
uuid = "996a588d-648d-4e1f-a8f0-a84b347e47b1"
version = "1.5.0"

    [deps.LogDensityProblemsAD.extensions]
    LogDensityProblemsADEnzymeExt = "Enzyme"
    LogDensityProblemsADFiniteDifferencesExt = "FiniteDifferences"
    LogDensityProblemsADForwardDiffBenchmarkToolsExt = ["BenchmarkTools", "ForwardDiff"]
    LogDensityProblemsADForwardDiffExt = "ForwardDiff"
    LogDensityProblemsADReverseDiffExt = "ReverseDiff"
    LogDensityProblemsADTrackerExt = "Tracker"
    LogDensityProblemsADZygoteExt = "Zygote"

    [deps.LogDensityProblemsAD.weakdeps]
    BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"
weakdeps = ["ChainRulesCore", "ChangesOfVariables", "InverseFunctions"]

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "0d097476b6c381ab7906460ef1ef1638fbce1d91"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.2"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "c88a4afe1703d731b1c4fdf4e3c7e77e3b176ea2"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.165"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Lux]]
deps = ["ADTypes", "Adapt", "ChainRulesCore", "ConcreteStructs", "Functors", "LinearAlgebra", "LuxCore", "LuxDeviceUtils", "LuxLib", "MacroTools", "Markdown", "Optimisers", "PackageExtensionCompat", "Random", "Reexport", "Setfield", "SparseArrays", "Statistics", "TruncatedStacktraces", "WeightInitializers"]
git-tree-sha1 = "42569bb572177026197dec4eca9834bc2fb94d7e"
uuid = "b2108857-7c20-44ae-9111-449ecde12c47"
version = "0.5.6"

    [deps.Lux.extensions]
    LuxComponentArraysExt = "ComponentArrays"
    LuxComponentArraysReverseDiffExt = ["ComponentArrays", "ReverseDiff"]
    LuxFluxTransformExt = "Flux"
    LuxLuxAMDGPUExt = "LuxAMDGPU"
    LuxLuxCUDAExt = "LuxCUDA"
    LuxTrackerExt = "Tracker"
    LuxZygoteExt = "Zygote"

    [deps.Lux.weakdeps]
    ComponentArrays = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
    FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
    Flux = "587475ba-b771-5e3f-ad9e-33799f191a9c"
    LuxAMDGPU = "83120cb1-ca15-4f04-bf3b-6967d2e6b60b"
    LuxCUDA = "d0bbae9a-e099-4d5b-a835-1c6931763bda"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.LuxCore]]
deps = ["Functors", "Random", "Setfield"]
git-tree-sha1 = "6feb02e23f6d70f407af97ca270b3a337af5ae0f"
uuid = "bb33d45b-7691-41d6-9220-0943567d0623"
version = "0.1.6"

[[deps.LuxDeviceUtils]]
deps = ["Adapt", "ChainRulesCore", "Functors", "LuxCore", "PackageExtensionCompat", "Preferences", "Random", "SparseArrays"]
git-tree-sha1 = "6adf6442934e8cff3f1d1af4b0e2fcd3358b6d72"
uuid = "34f89e08-e1d5-43b4-8944-0b49ac560553"
version = "0.1.8"

    [deps.LuxDeviceUtils.extensions]
    LuxDeviceUtilsFillArraysExt = "FillArrays"
    LuxDeviceUtilsLuxAMDGPUExt = "LuxAMDGPU"
    LuxDeviceUtilsLuxCUDAExt = "LuxCUDA"
    LuxDeviceUtilsMetalExt = "Metal"
    LuxDeviceUtilsZygoteExt = "Zygote"

    [deps.LuxDeviceUtils.weakdeps]
    FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
    LuxAMDGPU = "83120cb1-ca15-4f04-bf3b-6967d2e6b60b"
    LuxCUDA = "d0bbae9a-e099-4d5b-a835-1c6931763bda"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.LuxLib]]
deps = ["ChainRulesCore", "KernelAbstractions", "Markdown", "NNlib", "PackageExtensionCompat", "PrecompileTools", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "39e07c7373cb1dbff67591c5c2dd30c5a6cd7560"
uuid = "82251201-b29d-42c6-8e01-566dec8acb11"
version = "0.3.6"

    [deps.LuxLib.extensions]
    LuxLibForwardDiffExt = "ForwardDiff"
    LuxLibLuxCUDAExt = "LuxCUDA"
    LuxLibLuxCUDATrackerExt = ["LuxCUDA", "Tracker"]
    LuxLibReverseDiffExt = "ReverseDiff"
    LuxLibTrackerExt = "Tracker"

    [deps.LuxLib.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LuxCUDA = "d0bbae9a-e099-4d5b-a835-1c6931763bda"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.MCMCChains]]
deps = ["AbstractMCMC", "AxisArrays", "Dates", "Distributions", "Formatting", "IteratorInterfaceExtensions", "KernelDensity", "LinearAlgebra", "MCMCDiagnosticTools", "MLJModelInterface", "NaturalSort", "OrderedCollections", "PrettyTables", "Random", "RecipesBase", "Statistics", "StatsBase", "StatsFuns", "TableTraits", "Tables"]
git-tree-sha1 = "8778ea7283a0bf0d7e507a0235adfff38071493b"
uuid = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
version = "6.0.3"

[[deps.MCMCDiagnosticTools]]
deps = ["AbstractFFTs", "DataAPI", "DataStructures", "Distributions", "LinearAlgebra", "MLJModelInterface", "Random", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "3e6db72c2ab9cadfa3278ff388473a01fc0cfb9d"
uuid = "be115224-59cd-429b-ad48-344e309966f0"
version = "0.3.5"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "eb006abbd7041c28e0d16260e50a24f8f9104913"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2023.2.0+0"

[[deps.MLJModelInterface]]
deps = ["Random", "ScientificTypesBase", "StatisticalTraits"]
git-tree-sha1 = "03ae109be87f460fe3c96b8a0dbbf9c7bf840bd5"
uuid = "e80e1ace-859a-464e-9ed9-23947d8ae3ea"
version = "1.9.2"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MLUtils]]
deps = ["ChainRulesCore", "Compat", "DataAPI", "DelimitedFiles", "FLoops", "NNlib", "Random", "ShowCases", "SimpleTraits", "Statistics", "StatsBase", "Tables", "Transducers"]
git-tree-sha1 = "3504cdb8c2bc05bde4d4b09a81b01df88fcbbba0"
uuid = "f1d291b0-491e-4a28-83b9-f70985020b54"
version = "0.4.3"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathML]]
deps = ["AbstractTrees", "DocStringExtensions", "EzXML", "IfElse", "SpecialFunctions", "Statistics", "Symbolics"]
git-tree-sha1 = "12b0736d966c45d8621e895671266f05c881dda1"
uuid = "abcecc63-2b08-419c-80c4-c63dca6fa478"
version = "0.1.15"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "78f6e33434939b0ac9ba1df81e6d005ee85a7396"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "2.1.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "629afd7d10dbc6935ec59b32daeb33bc4460a42e"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.4"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.ModelingToolkit]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "Compat", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffRules", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "ForwardDiff", "FunctionWrappersWrappers", "Graphs", "IfElse", "InteractiveUtils", "JuliaFormatter", "JumpProcesses", "LabelledArrays", "Latexify", "Libdl", "LinearAlgebra", "MLStyle", "MacroTools", "NaNMath", "OrdinaryDiffEq", "PrecompileTools", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Serialization", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils", "Symbolics", "URIs", "UnPack", "Unitful"]
git-tree-sha1 = "87a45c453295c1640a1fd011a85e42fda7a78d15"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "8.70.0"

    [deps.ModelingToolkit.extensions]
    MTKDeepDiffsExt = "DeepDiffs"

    [deps.ModelingToolkit.weakdeps]
    DeepDiffs = "ab62b9b5-e342-54a8-a765-a90f495de1a6"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultiKDE]]
deps = ["Distributions"]
git-tree-sha1 = "8a48220d6d724d2b2bf3aaa95384eade42905031"
uuid = "77543b7f-bd95-4024-91c1-46775346e0e7"
version = "0.1.1"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "6c2e970692b6f4fed2508865c43a0f67f3820cf4"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.2"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "68bf5103e002c44adfd71fea6bd770b3f0586843"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.2"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "6985021d02ab8c509c841bb8b2becd3145a7b490"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.3.3"

[[deps.MvNormalCDF]]
deps = ["Distributions", "FillArrays", "LinearAlgebra", "Primes", "Random", "StatsBase"]
git-tree-sha1 = "e472ed20dfe9f68d9aebc02258dfe31d69deb8eb"
uuid = "37188c8d-bc69-4638-b057-733e744175ec"
version = "0.2.6"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NNlib]]
deps = ["Adapt", "Atomix", "ChainRulesCore", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "Pkg", "Random", "Requires", "Statistics"]
git-tree-sha1 = "3bc568de99214f72a76c7773ade218819afcc36e"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.9.7"

    [deps.NNlib.extensions]
    NNlibAMDGPUExt = "AMDGPU"
    NNlibCUDACUDNNExt = ["CUDA", "cuDNN"]
    NNlibCUDAExt = "CUDA"
    NNlibEnzymeCoreExt = "EnzymeCore"

    [deps.NNlib.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    cuDNN = "02a925ec-e4fe-4b08-9a7e-0d78e3d38ccd"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NameResolution]]
deps = ["PrettyPrint"]
git-tree-sha1 = "1a0fa0e9613f46c9b8c11eee38ebb4f590013c5e"
uuid = "71a1bf82-56d0-4bbc-8a3c-48b961074391"
version = "0.1.5"

[[deps.NamedArrays]]
deps = ["Combinatorics", "DataStructures", "DelimitedFiles", "InvertedIndices", "LinearAlgebra", "Random", "Requires", "SparseArrays", "Statistics"]
git-tree-sha1 = "6d42eca6c3a27dc79172d6d947ead136d88751bb"
uuid = "86f7a689-2022-50b4-a561-43c23ac3c673"
version = "0.10.0"

[[deps.NamedTupleTools]]
git-tree-sha1 = "90914795fc59df44120fe3fff6742bb0d7adb1d0"
uuid = "d9ec5142-1e00-5aa0-9d6a-321866360f50"
version = "0.14.3"

[[deps.NaturalSort]]
git-tree-sha1 = "eda490d06b9f7c00752ee81cfa451efe55521e21"
uuid = "c020b1a1-e9b0-503a-9c33-f039bfc54a85"
version = "1.0.0"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "2c3726ceb3388917602169bed973dbc97f1b51a8"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.13"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NoiseRobustDifferentiation]]
deps = ["IterativeSolvers", "LinearAlgebra", "LinearMaps", "Preconditioners", "SparseArrays"]
git-tree-sha1 = "108fc03c3419164898540545988d0df13fa5239b"
uuid = "470638dc-0858-4731-a73a-678bdc45695b"
version = "0.2.4"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "EnumX", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "PrecompileTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "e10debcea868cd6e51249e8eeaf191c25f68a640"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "1.10.1"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a12e56c72edee3ce6b96667745e6cbbe5498f200"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.23+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "963b004d15216f8129f6c0f7d187efa136570be0"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.7"

[[deps.Optimisers]]
deps = ["ChainRulesCore", "Functors", "LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "34205b1204cc83c43cd9cfe53ffbd3b310f6e8c5"
uuid = "3bd65402-5787-11e9-1adc-39752487f4e2"
version = "0.3.1"

[[deps.Optimization]]
deps = ["ADTypes", "ArrayInterface", "ConsoleProgressMonitor", "DocStringExtensions", "LinearAlgebra", "Logging", "LoggingExtras", "Pkg", "Printf", "ProgressLogging", "Reexport", "Requires", "SciMLBase", "SparseArrays", "TerminalLoggers"]
git-tree-sha1 = "f59193f79f8310e5eaad309ccfa8f2fa06b6703b"
uuid = "7f7a1694-90dd-40f0-9382-eb1efda571ba"
version = "3.19.1"

    [deps.Optimization.extensions]
    OptimizationEnzymeExt = "Enzyme"
    OptimizationFiniteDiffExt = "FiniteDiff"
    OptimizationForwardDiffExt = "ForwardDiff"
    OptimizationMTKExt = "ModelingToolkit"
    OptimizationReverseDiffExt = "ReverseDiff"
    OptimizationSparseDiffExt = ["SparseDiffTools", "Symbolics", "ReverseDiff"]
    OptimizationTrackerExt = "Tracker"
    OptimizationZygoteExt = "Zygote"

    [deps.Optimization.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseDiffTools = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.OptimizationBBO]]
deps = ["BlackBoxOptim", "Optimization", "Reexport"]
git-tree-sha1 = "939b19f2c213901aa544d1fd50e5f0065e66127c"
uuid = "3e6eede4-6085-4f62-9a71-46d9bc1eb92b"
version = "0.1.5"

[[deps.OptimizationOptimJL]]
deps = ["Optim", "Optimization", "Reexport", "SparseArrays"]
git-tree-sha1 = "76a57a8484370bb26ae82a0a330abfff96d561aa"
uuid = "36348300-93cb-4f02-beb5-3c3902f8871e"
version = "0.1.9"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLNLSolve", "SciMLOperators", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "16d4feaf7e17ac63cb6bb742b0a4dfd50d1ef458"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.57.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "b7c4f29f93b548caa58f703580f4d79ab753c8ac"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.21"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.PartialFunctions]]
git-tree-sha1 = "b3901ea034cfd8aae57a2fa0dde0b0ea18bad1cb"
uuid = "570af359-4316-4cb7-8c74-252c00c2016b"
version = "1.1.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "c7dc9720390fcc296bf757b3f833f9e41c68a086"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.7"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "ea78a2764f31715093de7ab495e12c0187f231d1"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.4"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preconditioners]]
deps = ["AlgebraicMultigrid", "LimitedLDLFactorizations", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "a02fc23058e6264b0f73523e99e2dd1a20232dff"
uuid = "af69fa37-3177-5a40-98ee-561f696e4fcd"
version = "0.6.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.PrettyPrint]]
git-tree-sha1 = "632eb4abab3449ab30c5e1afaa874f0b98b586e4"
uuid = "8162dcfd-2161-5ef2-ae6c-7681170c5f98"
version = "0.2.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "ee094908d720185ddbdc58dbe0c1cbe35453ec7a"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.7"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "4c9f306e5d6603ae203c2000dd460d81a5251489"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "00099623ffee15972c16111bcf84c58a0051257c"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.9.0"

[[deps.PumasQSP]]
deps = ["DataFrames", "DataSets", "DiffEqCallbacks", "Distributions", "JuliaSimModelOptimizer", "LinearAlgebra", "ModelingToolkit", "RecipesBase", "Reexport", "SBMLToolkit", "SciMLBase", "SciMLNLSolve", "SteadyStateDiffEq", "SymbolicUtils", "Symbolics", "UUIDs"]
git-tree-sha1 = "b1a84cfaf9489af39631c53a182f99efa0098309"
uuid = "5526fe93-3145-4cc2-889a-1fb5ef00f63e"
version = "2.4.5"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "7c29f0e8c575428bd84dc3c72ece5178caa67336"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.2+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9ebcd48c498668c7fa0e97a9cae873fbee7bfee1"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.1"

[[deps.QuasiMonteCarlo]]
deps = ["Distributions", "LatinHypercubeSampling", "LatticeRules", "LinearAlgebra", "Primes", "Random", "Sobol", "StatsBase"]
git-tree-sha1 = "b10f1a5345f14f6431712661925ae2d8c02780ee"
uuid = "8a4e6c94-4038-4cdc-81c3-7e6ffdb2a71b"
version = "0.2.19"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "552f30e847641591ba3f39fd1bed559b9deb0ef3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.6.1"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "d7087c013e8a496ff396bae843b1e16d9a30ede8"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.10"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "2b6d4a40339aa02655b1743f4cd7c03109f520c1"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.20"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.ReplMaker]]
deps = ["REPL", "Unicode"]
git-tree-sha1 = "f8bb680b97ee232c4c6591e213adc9c1e4ba0349"
uuid = "b873ce64-0db9-51f5-a568-4457d8e49576"
version = "0.2.7"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.ResourceContexts]]
deps = ["Logging"]
git-tree-sha1 = "0e9863272a09aff6579a987dd7fe3f94e32f6673"
uuid = "8d208092-d35c-4dd3-a0d7-8325f9cce6b4"
version = "0.2.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "06b5ac80ff1b88bd82df92c1c1875eea3954cd6e"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.20"

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[deps.SBML]]
deps = ["DocStringExtensions", "IfElse", "Libdl", "SBML_jll", "SparseArrays", "Unitful"]
git-tree-sha1 = "28d25ec1f0bdadd291056e4d421d363f4c159309"
uuid = "e5567a89-2604-4b09-9718-f5f78e97c3bb"
version = "1.5.0"

[[deps.SBMLToolkit]]
deps = ["Catalyst", "MathML", "SBML", "SymbolicUtils"]
git-tree-sha1 = "4f39460a922221c5864762ac85a569de3bce1f99"
uuid = "86080e66-c8ac-44c2-a1a0-9adaadfe4a4e"
version = "0.1.25"

[[deps.SBML_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "7567084bd3defc697c168b8b0ab3d28a9bc47648"
uuid = "bb12108a-f4ef-5f88-8ef3-0b33ff7017f1"
version = "5.20.1+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "0e270732477b9e551d884e6b07e23bb2ec947790"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.5"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "4b8586aece42bee682399c4c4aee95446aa5cd19"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.39"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "ChainRulesCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FillArrays", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "916b8a94c0d61fa5f7c5295649d3746afb866aff"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.98.1"

    [deps.SciMLBase.extensions]
    ZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLNLSolve]]
deps = ["DiffEqBase", "LineSearches", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "765b788339abd7d983618c09cfc0192e2b6b15fd"
uuid = "e9a6253c-8580-4d32-9898-8661bb511710"
version = "0.1.9"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "65c2e6ced6f62ea796af251eb292a0e131a3613b"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.6"

[[deps.ScientificTypesBase]]
git-tree-sha1 = "a8e18eb383b5ecf1b5e6fc237eb39255044fd92b"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "3.0.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "04bdff0b09c65ff3e06a05e3eb7b120223da3d39"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SetRounding]]
git-tree-sha1 = "d7a25e439d07a17b7cdf97eecee504c50fedf5f6"
uuid = "3cc68bcd-71a2-5612-b932-767ffbe40ab0"
version = "0.2.1"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.ShowCases]]
git-tree-sha1 = "7f534ad62ab2bd48591bdeac81994ea8c445e4a5"
uuid = "605ecd9f-84a6-4c9e-81e2-4798472b76a3"
version = "0.1.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "PackageExtensionCompat", "PrecompileTools", "Reexport", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "4d53b83af904049c493daaf2a225bcae994a3c59"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "0.1.20"
weakdeps = ["NNlib"]

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveNNlibExt = "NNlib"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sobol]]
deps = ["DelimitedFiles", "Random"]
git-tree-sha1 = "5a74ac22a9daef23705f010f72c81d6925b19df8"
uuid = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
version = "1.5.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "Tricks", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "42d131931906bf4f0af97a7113c8456d0a8aff9d"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.6.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SparseInverseSubset]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "91402087fd5d13b2d97e3ef29bbdf9d7859e678a"
uuid = "dc90abb0-5640-4711-901d-7e5b23a2fada"
version = "0.1.1"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpatialIndexing]]
git-tree-sha1 = "bacf5065cd7c0d6449b8bba6fa8e75b3087356b0"
uuid = "d4ead438-fe20-5cc5-a293-4fd39a41b74c"
version = "0.1.5"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StableRNGs]]
deps = ["Random", "Test"]
git-tree-sha1 = "3be7d49667040add7ee151fefaf1f8c04c8c8276"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.0"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "f295e0a1da4ca425659c57441bcb59abb035a4bc"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.8"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "03fec6800a986d191f64f5c0996b59ed526eda25"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.4.1"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "0adf069a2a490c47273727e029371b31d44b72b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.5"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "30b9236691858e13f167ce829490a68e1a597782"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "3.2.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "9115a29e6c2cf66cf213ccc17ffd61e27e743b24"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.6"

[[deps.SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "2ca69f4be3294e4cd987d83d6019037d420d9fc1"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.16.1"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "b341540a647b39728b6d64eaeda82178e848f76e"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.62.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "f02eb61eb5c97b48c153861c72fbbfdddc607e06"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.17"

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "b765e46ba27ecf6b44faf70df40c57aa3a547dcb"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.StructArrays]]
deps = ["Adapt", "ConstructionBase", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "0a3db38e4cce3c54fe7a71f831cd7b6194a54213"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.16"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "b208d5379feee6c1c9f9660e4294e77e40946884"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.20.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "04777432d74ec5bc91ca047c9e0e0fd7f81acdb6"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "2f3fa844bcd33e40d8c29de5ee8dded7a0a70422"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.4.0"

[[deps.Symbolics]]
deps = ["ArrayInterface", "Bijections", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "d8b3a50cf6ccb19a090c0bb6d89fd4f1b576fcfb"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.6.0"

    [deps.Symbolics.extensions]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a1f34829d5ac0ef499f6d84428bd6b4c71f02ead"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaylorSeries]]
deps = ["IntervalArithmetic", "LinearAlgebra", "Markdown", "Requires", "SparseArrays"]
git-tree-sha1 = "b123e0614b0f366e2073c4da7b1367cebda3f579"
uuid = "6aa5eb33-94cf-58f4-a9d0-e4b2c4fc25ea"
version = "0.14.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TerminalLoggers]]
deps = ["LeftChildRightSiblingTrees", "Logging", "Markdown", "Printf", "ProgressLogging", "UUIDs"]
git-tree-sha1 = "f133fab380933d042f6796eda4e130272ba520ca"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.7"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "34e6bcf36b9ed5d56489600cf9f3c16843fa2aa2"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.11"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.ToeplitzMatrices]]
deps = ["AbstractFFTs", "DSP", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "a097662c0fd28143bb714e4527ac2dcb7520212e"
uuid = "c751599d-da0a-543b-9d20-d0a503d91d24"
version = "0.8.2"
weakdeps = ["StatsBase"]

    [deps.ToeplitzMatrices.extensions]
    ToeplitzMatricesStatsBaseExt = "StatsBase"

[[deps.Tokenize]]
git-tree-sha1 = "90538bf898832b6ebd900fa40f223e695970e3a5"
uuid = "0796e94c-ce3b-5d07-9a54-7f471281c624"
version = "0.5.25"

[[deps.Tracker]]
deps = ["Adapt", "DiffRules", "ForwardDiff", "Functors", "LinearAlgebra", "LogExpFunctions", "MacroTools", "NNlib", "NaNMath", "Optimisers", "Printf", "Random", "Requires", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "994c4977c13e7e2525566e193f82374021fab5dd"
uuid = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
version = "0.2.27"
weakdeps = ["PDMats"]

    [deps.Tracker.extensions]
    TrackerPDMatsExt = "PDMats"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "53bd5978b182fa7c57577bdb452c35e5b4fb73a5"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.78"

    [deps.Transducers.extensions]
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Trapz]]
git-tree-sha1 = "79eb0ed763084a3e7de81fe1838379ac6a23b6a0"
uuid = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
version = "2.0.3"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "31eedbc0b6d07c08a700e26d31298ac27ef330eb"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.19"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.Turing]]
deps = ["AbstractMCMC", "AdvancedHMC", "AdvancedMH", "AdvancedPS", "AdvancedVI", "BangBang", "Bijectors", "DataStructures", "Distributions", "DistributionsAD", "DocStringExtensions", "DynamicPPL", "EllipticalSliceSampling", "ForwardDiff", "Libtask", "LinearAlgebra", "LogDensityProblems", "LogDensityProblemsAD", "MCMCChains", "NamedArrays", "Printf", "Random", "Reexport", "Requires", "SciMLBase", "Setfield", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Tracker"]
git-tree-sha1 = "5e9b655060490d1a16cbe569acb276578ace48aa"
uuid = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"
version = "0.28.3"

    [deps.Turing.extensions]
    TuringDynamicHMCExt = "DynamicHMC"
    TuringOptimExt = "Optim"

    [deps.Turing.weakdeps]
    DynamicHMC = "bbc10e6e-7c05-544b-b16e-64fede858acb"
    Optim = "429524aa-4258-5aef-a3af-852621145aeb"

[[deps.URIs]]
git-tree-sha1 = "b7a5e99f24892b6824a954199a45e9ffcc1c70f0"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a72d22c7e13fe2de562feda8645aa134712a87ee"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.17.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "21c8fc7cd598ef49f11bc9e94871f5d7740e34b9"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.5"

[[deps.UnsafeAtomics]]
git-tree-sha1 = "6331ac3440856ea1988316b46045303bef658278"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.2.1"

[[deps.UnsafeAtomicsLLVM]]
deps = ["LLVM", "UnsafeAtomics"]
git-tree-sha1 = "323e3d0acf5e78a56dfae7bd8928c989b4f3083e"
uuid = "d80eeb9a-aca5-4d75-85e5-170c8b632249"
version = "0.1.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "b182207d4af54ac64cbc71797765068fdeff475d"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.64"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WeightInitializers]]
deps = ["PartialFunctions", "Random", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "214ad59a8ece2118a1e70a07a32cef784a5c0704"
uuid = "d49dbf32-c5c2-4618-8acc-27bb2598ef2d"
version = "0.1.1"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "24b81b59bd35b3c42ab84fa589086e19be919916"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.11.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cf2c7de82431ca6f39250d2fc4aacd0daa1675c0"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.4+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "e6330e4b731a6af7959673621e91645eb1356884"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.9"

[[deps.ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "f492b7fe1698e623024e873244f10d89c95c340a"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.10.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "977aed5d006b840e2e40c0b48984f7463109046d"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.3"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─41f23ebc-cc05-46bc-9879-b986bb16dbc1
# ╟─a1e18f74-7f8f-4a05-b7d7-22b914fda4b8
# ╟─bbbfb448-771b-4b0c-b7e9-627b2fb2d2b7
# ╟─f960d305-eef9-46e8-93d5-b7980f70a159
# ╠═2bd1af72-2f3e-471b-bb28-32e38c1a32bd
# ╟─79d9a06f-40ab-44c8-b701-7d877e70eb81
# ╟─56c50a67-c7f8-4009-9f4d-e624097b1875
# ╠═30f6f698-7cbe-4f4c-8791-6c642c4fd816
# ╟─a96e0591-0cca-4695-9a43-5f7e5204e83c
# ╠═8ebf0e5e-9980-4511-805f-b99d14803077
# ╟─38681cc1-e4b1-42c4-8836-79c9f145029b
# ╟─2266bace-b5a9-40f0-85df-ed37ee4d5dd5
# ╟─40611a30-696e-415b-b082-1ac96eeefd0d
# ╠═59c3d5bf-8dc2-4f37-b7c0-1eb4d8f339ba
# ╟─995e26f2-4aee-4239-a697-3ee895301bc8
# ╠═a3d1b87f-2246-4212-ab2d-4370e9e5278f
# ╟─81beeae7-7599-45fb-877f-baf2c09c16fb
# ╠═617615bd-6aa1-4994-815d-65f68d0ac474
# ╟─721fb15d-853b-4ee1-89ac-154b28590f53
# ╠═6eb29af1-29e2-4af6-baad-08d745b8f7c5
# ╠═0d3aae0a-0c7e-4dd9-bee0-7496a54a96d6
# ╠═f6e787a8-fb1d-4802-812b-966b1a4bd334
# ╟─6faa7f93-144b-4d1f-9aef-87679cb5025a
# ╠═c42e63de-ccec-4434-b541-23ae3b378e72
# ╟─af6bf403-35a3-47dc-9312-e9c67ce6c20d
# ╟─a94c1bd6-2fbc-4412-b36d-7acc75864335
# ╟─09663cb0-8f4d-43cd-9d18-2a48b04f782e
# ╟─904c373d-a38b-4130-8b1c-dc044631b72d
# ╟─7ca4fe49-3f21-4c30-91d4-459ff7642368
# ╟─7ce39e9f-6fd7-4987-94f7-6c1da9622183
# ╟─0d95c7be-7679-4435-ab8a-0f69c59fd747
# ╟─937e5c93-ffd1-4b67-8350-37d83236d082
# ╟─c230fc3c-17e4-4426-abfe-139311d4a26e
# ╟─3faede93-38f6-432b-b39c-553bc89260ab
# ╟─a2f66f61-bfab-4879-920c-ec120a226be7
# ╟─3f609536-e997-4dd0-9843-db7946e00feb
# ╟─f1b7b3c6-bdae-4314-8f1e-fac63d908c6a
# ╟─1d263fee-1eda-47f5-a0c0-a02e956c0a09
# ╟─6752199d-2fc3-4e4e-81ba-0788ffb4a7f9
# ╟─37f4c431-651d-4028-95c2-3a304ecca45e
# ╟─12954f4f-ccea-4323-8b70-47ef378b4064
# ╟─b573d05a-b265-4f7c-b434-9a7d22c68464
# ╟─51194b3c-72ea-418f-9f06-66531cb50960
# ╟─be61f182-d02d-41b0-8648-ee896d819ce1
# ╟─70819d97-7187-4b92-87b7-29ae202fa0f6
# ╟─07d35ca9-40b3-43f5-8bd1-eb4c06710a3e
# ╟─79ecafd1-745c-4889-a1b7-d9eaceb9a5d7
# ╟─88592d7f-8d41-4943-b381-36f9aedb9b8f
# ╟─7f5842a5-c184-4f59-92cd-fb616e546462
# ╟─7d9469a1-cd1d-4696-9661-f3ad53cb72fa
# ╟─4d8631a5-b770-4659-b7f5-7d37511d730f
# ╟─1362250a-c20a-4709-9c2c-55bdbf8b3e9d
# ╟─398ad035-21c7-4055-a206-806c69bcaa1a
# ╟─d91ec775-e2ef-4b78-96e3-eeb789ca7f55
# ╟─d4837b7c-2d4f-4289-9f06-43922856fb1d
# ╟─34d5c411-ce58-4d7c-9632-7f4b91a1d480
# ╟─a0b29e6b-5f52-4e75-85bb-cf3a6e0b2d20
# ╟─60a1677c-ffb2-4b73-97a0-a3cd2c2a5580
# ╟─f580f2f6-60a0-46e3-b8ff-57090a607a03
# ╟─6ba09f38-313d-4e90-a83b-f902f42a4ea4
# ╟─265dff45-b8e6-4e28-aac7-b557ff52eddb
# ╟─75fbd911-9995-48ce-b61c-2f9083b5aa16
# ╟─0b376a7a-7906-438c-9738-56390a03af54
# ╟─eeb8a453-e9af-400e-a7c1-5c994a97b2f8
# ╠═eea3747a-23b9-4653-a789-2ab7f6f3055f
# ╟─b74abd28-7fd0-4e01-b821-c53acef1f715
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
