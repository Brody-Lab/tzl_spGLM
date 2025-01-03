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
	Options(csvpath, datapath)
"""
function Options(csvpath::String, datapath::String)
	dict = dictionary(csvpath,1)
	dict["datapath"] = datapath
	Options(dict)
end

"""
	dictionary(csvpath, row)

RETURN a `Dict` containing hyperparameters

ARGUMENT
-`csvpath`: the absolute path to a comma-separated values (CSV) file
-`row`: the row of the CSV to be considered
"""
function dictionary(csvpath::String, row::Integer)
	df = DataFrames.DataFrame(CSV.File(csvpath))
	dfrow = df[row,:]
	Dict((name=>dfrow[name] for name in names(dfrow))...)
end

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
						if ismissing(options[String(fieldname)])
							""
						else
							convert(typeof(defaultvalue), options[String(fieldname)])
						end
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
	if options.outputpath != "placeholder"
		!isdir(options.outputpath) && mkpath(options.outputpath)
		@assert isdir(options.outputpath)
	end
	return options
end

"""
RETURN a `Dict` containing the default values of the fixed hyperparameters ("options") of the model 
"""
dictionary_default_options() = SPGLM.dictionary(SPGLM.Options())

"""
	namedoptions(optionsname)

Load the a set of hyperparameters saved in `/options.csv`
"""
function dictionary_named_options(optionsname::String)
	csvpath = joinpath(dirname(@__DIR__),"options.csv")
	df = DataFrames.DataFrame(CSV.File(csvpath))
	row = df[df.name .== optionsname,:]
	select!(row, Not(:name))
	dict = Dict(pairs(eachcol(row)))
	Dict(string(k) => v[1] for (k, v) in dict)
end

"""
	Options(options, fieldname, value)

Replace the value of a field in an instance of the struct `Options`

ARGUMENT
-`options`: a structure containing the fixed hyperparameters of the model
-`fieldname`: the field to be replaced
-`value`: the new value of the field
"""
function Options(options::Options, fieldname::Symbol, value)
	Options(Dict((String(name)=>(name==fieldname ? value : getfield(options,name)) for name in fieldnames(SPGLM.Options))...))
end

"""
	loadtrials(options)

RETURN a vector of objects of the composite type `trials`

ARGUMENT
-`options`: a structure containing fixed hyperparameters
"""
function loadtrials(options::Options)
	if (options.projectname == "") || (options.projectname == "manuscript2023a")
		loadtrials_manuscript2023a(options)
	elseif options.projectname == "auditory"
		loadtrials_auditory(options)
	end
end

"""
	loadtrials_manuscript2023a(options)

Load trials for the project "manuscript2023a"

RETURN a vector of objects of the composite type `Trial`

ARGUMENT
-`options`: a structure containing fixed hyperparameters
"""
function loadtrials_manuscript2023a(options::Options)
	file = matopen(options.datapath)
	Cell = read(file, "Cell")
	if haskey(Cell, "Trials_path")
		Trialsfile = matopen(Cell["Trials_path"])
		Trials = read(Trialsfile)["Trials"]
		close(Trialsfile)
	else
		Trials = read(file, "Trials")
	end
	close(file)
	if !isempty(options.trial_indices_path)
		if options.trial_indices_path == "do_not_exclude"
			trialindices = collect(1:length(Trials["pokedR"]))
		else
			trialindices = DataFrames.DataFrame(CSV.File(options.trial_indices_path; header=0))[:,1]
		end
	else
		trialindices = (.!Trials["violated"]) .&
					   	(Trials["trial_type"] .== "a") .&
					   	.!isnan.(Trials["pokedR"]) .&
					   	.!isnan.(Trials["stateTimes"]["cpoke_in"]) .&
						.!isnan.(Trials["stateTimes"]["spoke"]) .&
						(Trials["stateTimes"]["cpoke_out"] .> Trials["stateTimes"]["clicks_on"])
		if haskey(Trials, "responded")
			trialindices = trialindices .& Trials["responded"]
		end
		trialindices = findall(vec(trialindices))
	end
	Na = floor(Int, options.time_in_trial_begin_s/options.dt)
	Npre = ceil(Int, (options.bfs_postspike_end_s-options.bfs_postspike_begin_s)/options.dt)
	pre_binedges = -Npre*options.dt:options.dt:0
	spiketimes_s = vec(Cell["spiketimes_s"])
	if options.reference_event == "stereoclick"
		first_reference_time_s = Trials["leftBups"][trialindices[1]][1] + Trials["stateTimes"]["clicks_on"][trialindices[1]]
		last_reference_time_s = Trials["leftBups"][trialindices[end]][1] + Trials["stateTimes"]["clicks_on"][trialindices[end]]
	else
		first_reference_time_s = Trials["stateTimes"][options.reference_event][trialindices[1]]
		last_reference_time_s = Trials["stateTimes"][options.reference_event][trialindices[end]]
	end
	usepose = (options.pose_filepath != "")
	if usepose
		file = matopen(options.pose_filepath)
		pose = vec(vec.(read(file, "pose")))
		pose_time_s = vec(read(file, "pose_time_s"))
		close(file)
	end
	map(trialindices) do i
		if options.reference_event == "stereoclick"
			reference_time_s = Trials["leftBups"][i][1] + Trials["stateTimes"]["clicks_on"][i]
		else
			reference_time_s = Trials["stateTimes"][options.reference_event][i]
		end
		t₀ = reference_time_s + Na*options.dt
		if isnan(options.time_in_trial_end_s) && !isempty(options.trim_after_event) # the trial is truncated with respect to only `trim_after_event` and not `reference_event`
			Nb = floor(Int, (Trials["stateTimes"][options.trim_after_event][i]-t₀)/options.dt)
		elseif !isnan(options.time_in_trial_end_s) && !isempty(options.trim_after_event) # the trial is truncated with respect to both `trim_after_event` and `reference_event`
			lasttimestep1 = floor(Int, (Trials["stateTimes"][options.trim_after_event][i]-t₀)/options.dt)
			lasttimestep2 = ceil(Int, options.time_in_trial_end_s/options.dt)
			Nb = min(lasttimestep1,lasttimestep2)
			@assert Nb > 0
		else # the trial is truncated with respect to only `reference_event` and not `trim_after_event`
			Nb = ceil(Int, options.time_in_trial_end_s/options.dt)
		end
		binedges_s = (Na*options.dt):options.dt:(Nb*options.dt)
		hist = StatsBase.fit(Histogram, spiketimes_s, (reference_time_s .+ binedges_s), closed=:left)
		y = hist.weights
		ypre = StatsBase.fit(Histogram, spiketimes_s, (reference_time_s .+ pre_binedges), closed=:left).weights
		if usepose # assumes sampling rate of the pose signals is lower than that of the neural
			ntimesteps = length(binedges_s)-1
			pose_sorted = collect(fill(NaN,ntimesteps) for pose in pose)
			for j = 1:ntimesteps
				index = findfirst(pose_time_s .>= (reference_time_s+binedges_s[j]))
				for k in eachindex(pose)
					pose_sorted[k][j] = sum(pose[k][index])
				end
			end
		else
			pose_sorted = collect(zeros(0) for i=1:length(trialindices))
		end
		@assert !isnan(Trials["stateTimes"]["cpoke_out"][i])
		fixation_timestep = ceil(Int, (Trials["stateTimes"]["cpoke_in"][i] - t₀)/options.dt)
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
				first_reference_time_s=first_reference_time_s,
				fixation_timestep = fixation_timestep,
				γ = Trials["gamma"][i],
				last_reference_time_s=last_reference_time_s,
				movement_timestep = movement_timestep,
				pose=pose_sorted,
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
Model(options::Options, trials::Vector{<:Trial}) = Model(options,trials, baselineweights(options,trials))

function Model(options::Options, trials::Vector{<:Trial}, 𝐰_baseline::Vector{<:AbstractFloat})
	𝐲 = vcat((trial.y for trial in trials)...)
	T = length(𝐲)
	𝐗 = Array{typeof(1.0)}(undef,T,0)
	if !isempty(𝐰_baseline)
		𝐛 = baselineinput(options,trials)*𝐰_baseline
		𝐗 = hcat(𝐗, vcat((fill(b,trial.T) for (b,trial) in zip(𝐛,trials))...)) # using `hcat` for type
	else
		𝐛 = collect(NaN for trial in trials)
	end
	setnames = SPGLM.basis_function_sets()
	basissets = collect(BasisFunctionSet(setname, options, trials) for setname in setnames)
	indices = UnitRange{Int}[]
	k = 0
	for inputname in fieldnames(SPGLM.WeightIndices)
		if inputname==:baseline
			i = Int(!isempty(𝐰_baseline))
			indices = vcat(indices, [k .+ (1:i)])
			k += i
		elseif getfield(options, Symbol("input_"*String(inputname)))
			if inputname == :pose
				𝐗add = vcat((hcat(trial.pose...) for trial in trials)...)
			else
				basisset = first(basissets[setnames .== SPGLM.match_input_to_basis(inputname)])
				𝐗add = SPGLM.inputs_each_timestep(basisset, inputname, trials)
			end
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
			baseline=𝐛,
			basissets=basissets,
			trials=trials,
			weightindices = weightindices,
			𝐰_baseline = 𝐰_baseline,
			𝐗=𝐗,
			𝐲=𝐲)
end
