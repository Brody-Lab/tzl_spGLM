"""
	Options(csvpath, row)

RETURN a struct containing the fixed hyperparameters of the model

ARGUMENT
-`csvpath`: the absolute path to a comma-separated values (CSV) file
-`row`: the row of the CSV to be considered
"""
Options(csvpath::String) = Options(csvpath::String, 1)
Options(csvpath::String, row::Integer) = Options(DataFrames.DataFrame(CSV.File(csvpath)), row)
Options(df::DataFrames.DataFrame, row::Integer) = Options(df[row,:])
Options(dfrow::DataFrames.DataFrameRow) = Options(Dict((name=>dfrow[name] for name in names(dfrow))...))

"""
	Options(options::Dict)

RETURN a struct containing the fixed hyperparameters of the model

ARGUMENT
-`options`: a dictionary
"""
function Options(options::Dict)
	keyset = keys(options)
	defaults = Options()
	pathfields = (:datapath, :outputpath)
	entries = 	map(fieldnames(Options)) do fieldname
					defaultvalue = getfield(defaults,fieldname)
					if String(fieldname) ∈ keyset
						convert(typeof(defaultvalue), options[String(fieldname)])
					elseif fieldname ∈ pathfields
						convert(typeof(defaultvalue), options[String(fieldname)])
					else
						defaultvalue
					end
				end
	options = Options(entries...)
	for fieldname in pathfields
		@assert !isempty(getfield(options, fieldname))
	end
	@assert isfile(options.datapath)
	!isdir(options.outputpath) && mkpath(options.outputpath)
	@assert isdir(options.outputpath)
	return options
end

"""
	loadtrials(options)

RETURN a vector of objects of the composite type `trials`

ARGUMENT
-`options`: a structure containing fixed hyperparameters
"""
function loadtrials(options::Options)
	matfile = read(matopen(options.datapath))
	Cell = matfile["Cell"]
	Trials = matfile["Trials"]
	trialindices = (.!Trials["violated"]) .& (Trials["trial_type"] .== "a") .& (Trials["stateTimes"]["cpoke_out"] .> Trials["stateTimes"]["clicks_on"])
	trialindices = vec(trialindices)
	Na = floor(Int, options.time_in_trial_begin_s/options.dt)
	Nb = ceil(Int, options.time_in_trial_end_s/options.dt)
	default_binedges_s = (Na*options.dt):options.dt:(Nb*options.dt)
	default_N = length(default_binedges_s)
	Npre = ceil(Int, (options.bfs_postspike_end_s-options.bfs_postspike_begin_s)/options.dt)
	pre_binedges = -Npre*options.dt:options.dt:0
	spiketimes_s = vec(Cell["spiketimes_s"])
	map(findall(trialindices)) do i
		if options.reference_event == "stereoclick"
			reference_time_s = Trials["leftBups"][i][1] + Trials["stateTimes"]["clicks_on"][i]
		else
			reference_time_s = Trials["stateTimes"][options.reference_event][i]
		end
		if !isempty(options.trim_after_event)
			lasttimestep = floor(Int, (Trials["stateTimes"][options.trim_after_event][i]-reference_time_s)/options.dt)
			if lasttimestep < default_N
				binedges_s = default_binedges_s[1:lasttimestep]
			end
		else
			binedges_s = default_binedges_s
		end
		hist = StatsBase.fit(Histogram, spiketimes_s, (reference_time_s .+ binedges_s), closed=:right)
		y = hist.weights
		ypre = StatsBase.fit(Histogram, spiketimes_s, (reference_time_s .+ pre_binedges), closed=:right).weights
		@assert !isnan(Trials["stateTimes"]["cpoke_out"][i])
		t₀ = reference_time_s + binedges_s[1]
		movement_timestep = ceil(Int, (Trials["stateTimes"]["cpoke_out"][i] - t₀)/options.dt)
		leftclicks_s = Trials["leftBups"][i] .+ Trials["stateTimes"]["clicks_on"][i] .- t₀
		rightclicks_s = Trials["rightBups"][i] .+ Trials["stateTimes"]["clicks_on"][i] .- t₀
		if typeof(leftclicks_s)<:AbstractFloat
			leftclicks_s = [leftclicks_s]
		end
		if typeof(rightclicks_s)<:AbstractFloat
			rightclicks_s = [rightclicks_s]
		end
		@assert leftclicks_s[1] == rightclicks_s[1]
		stereoclick_timestep = ceil(Int, leftclicks_s[1]/options.dt)
		leftclicks_timestep = ceil.(Int, leftclicks_s[2:end]./options.dt)
		rightclicks_timestep = ceil.(Int, rightclicks_s[2:end]./options.dt)
		response_timestep = ceil(Int, (Trials["stateTimes"]["spoke"][i] - t₀)/options.dt)
		clicks_timestep = [stereoclick_timestep; leftclicks_timestep; rightclicks_timestep]
		clicks_source = [2; zeros(Int, length(leftclicks_timestep)); ones(Int, length(rightclicks_timestep))]
		Trial(choice = Trials["pokedR"][i] == 1,
				clicks_source=clicks_source,
				clicks_timestep=clicks_timestep,
				γ = Trials["gamma"][i],
				movement_timestep = movement_timestep,
				reference_time_s=reference_time_s,
				response_timestep = response_timestep,
				timesteps_s=binedges_s[2:end],
				trialindex=i,
				y=y,
				ypre=ypre)
	end
end

"""
	Model(csvpath, row)

RETURN a struct containing data, parameters, and hyperparameters

ARGUMENT
-`csvpath`: the absolute path to a comma-separated values (CSV) file
-`row`: the row of the CSV to be considered
"""
Model(csvpath::String, row::Integer) = Model(Options(csvpath, row))
Model(csvpath::String) = Model(csvpath,1)
Model(options::Options) = Model(options, loadtrials(options))

"""
    Model(options, trialsets)

RETURN a struct containing data, parameters, and hyperparameters of a factorial hidden Markov drift-diffusion model

ARGUMENT
-`options`: a struct containing the fixed hyperparameters of the model
-`trialsets`: data used to constrain the model
"""
function Model(options::Options, trials::Vector{<:Trial})
	𝐲 = vcat((trial.y for trial in trials)...)
	T = length(𝐲)
	𝐗 = Array{typeof(1.0)}(undef,T,0)
	setnames = SPGLM.basis_function_sets()
	basissets = collect(SPGLM.BasisFunctionSet(setname, options) for setname in setnames)
	indices = UnitRange{Int}[]
	k = 0
	for inputname in fieldnames(SPGLM.WeightIndices)
		if getfield(options, Symbol("input_"*String(inputname)))
			basisset = first(basissets[setnames .== SPGLM.match_input_to_basis(inputname)])
			𝐗add = SPGLM.inputs_each_timestep(basisset, inputname, trials)
			𝐗 = hcat(𝐗, 𝐗add)
			N = size(𝐗add,2)
			indices = vcat(indices, [k .+ (1:N)])
			k += N
		else
			indices = vcat(indices, [1:0])
		end
	end
	weightindices = SPGLM.WeightIndices(indices...)
	SPGLM.Model(options=options,
			basissets=basissets,
			trials=trials,
			weightindices = weightindices,
			𝐗=𝐗,
			𝐲=𝐲)
end
