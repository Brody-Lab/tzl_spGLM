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
"""
function loadtrials(options::Options)
	matfile = read(matopen(options.datapath))
	Cell = matfile["Cell"]
	Trials = matfile["Trials"]
	trialindices = (.!Trials["violated"]) .& (Trials["trial_type"] .== "a") .& (Trials["stateTimes"]["cpoke_out"] .> Trials["stateTimes"]["clicks_on"])
	trialindices = vec(trialindices)
	ntimesteps = convert(Integer, (options.reference_end_s - options.reference_begin_s)/options.dt) # intentionally to trigger `InexactError` if `ntimesteps` is not an integer
	binedges_s = options.reference_begin_s:options.dt:(ntimesteps*options.dt)
	reference_times_s = vec(Trials["stateTimes"][options.reference_event][trialindices])
	spiketimes_s = vec(Cell["spiketimes_s"])
	map(findall(trialindices)) do i
		reference_time_s = Trials["stateTimes"][options.reference_event][i]
		hist = StatsBase.fit(Histogram, spiketimes_s, (reference_time_s .+ binedges_s), closed=:right)
		y = hist.weights
		@assert !isnan(Trials["stateTimes"]["cpoke_out"][i])
		if !isempty(Trials["leftBups"][i])
			@assert Trials["leftBups"][i][1] == Trials["rightBups"][i][1]
			stereoclick_time_s = Trials["leftBups"][i][1]+Trials["stateTimes"]["clicks_on"][i]-reference_time_s
		else
			stereoclick_time_s = NaN;
		end
		if typeof(Trials["leftBups"][i])<:AbstractFloat
			Lclick_times_s = ones(0)
		else
			Lclick_times_s = vec(Trials["leftBups"][i][2:end])
		end
		if typeof(Trials["rightBups"][i])<:AbstractFloat
			Rclick_times_s = ones(0)
		else
			Rclick_times_s = vec(Trials["rightBups"][i][2:end])
		end
		Trial(binedges_s=binedges_s,
			  choice = Trials["pokedR"][i] == 1,
			  γ = Trials["gamma"][i],
			  Lclick_times_s = Lclick_times_s .+ Trials["stateTimes"]["clicks_on"][i] .- reference_time_s,
			  movement_time_s = Trials["stateTimes"]["cpoke_out"][i] - reference_time_s,
			  stereoclick_time_s = stereoclick_time_s,
			  Rclick_times_s = Rclick_times_s .+ Trials["stateTimes"]["clicks_on"][i] .- reference_time_s,
			  reference_time_s=reference_time_s,
			  trialindex=i,
			  y=y)
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

"""
    Model(options, trialsets)

RETURN a struct containing data, parameters, and hyperparameters of a factorial hidden Markov drift-diffusion model

ARGUMENT
-`options`: a struct containing the fixed hyperparameters of the model
-`trialsets`: data used to constrain the model
"""
function Model(options::Options, trials::Vector{<:Trial})
	Φclick = temporal_basis_functions("click", options)
	Φpostspike = temporal_basis_functions("postspike", options)
	Φmovement = temporal_basis_functions("movement", options)
	Φfixation = temporal_basis_functions("fixation", options)
	Φdrift, 𝐗drift = drift_design_matrix(options, stereoclick_times_s, trialdurations, 𝐲)
	𝐲 = vcat((trial.spiketrains[n] for trial in trials)...)
	T = length(𝐲)
	empty = Array{typeof(1.0)}(undef,T,0)
	𝐗 = copy(empty)
	for fieldname in fieldnames(WeightIndices)
		indices = UnitRange{Int}[]
		k = 0
		if contains(String(fieldname), "click")
			if getfield(options, Symbol("include_"*String(fieldname)))
				𝐗add = click_inputs(Φclick, trials, laterality=fieldname)
				𝐗 = hcat(𝐗, 𝐗add)
				N = size(𝐗add,2)
				indices = vcat(indices, [k .+ (1:N)])
				k += N
			end
		end
		if contains(String(fieldname), "movement")
			if getfield(options, Symbol("include_"*String(fieldname)))
				𝐗add = movement_inputs(Φmovement, trials, laterality=fieldname)
				𝐗 = hcat(𝐗, 𝐗add)
				N = size(𝐗add,2)
				indices = vcat(indices, [k .+ (1:N)])
				k += N
			end
		end
	end
	𝐗stereoclick = options.include_stereoclick ? click_inputs(Φclick, trials, laterality=2) : empty
	𝐗leftclick = click_inputs(Φclick, trials, laterality= -1)
	𝐗rightclick = click_inputs(Φclick, trials, laterality= 1)
	𝐗leftchoice = movement_inputs(Φmovement, trials, laterality= -1)
	𝐗rightchoice = movement_inputs(Φmovement, trials, laterality= 1)
	𝐗fixation = fixation_inputs(Φfixation, trials)
	𝐗postspike = postspike_inputs(Φpostspike, trialdurations, 𝐲)
	𝐗=hcat(𝐗postspike, 𝐗drift, 𝐗fixation, 𝐗stereoclick, 𝐗leftclick, 𝐗rightclick, 𝐗leftchoice, 𝐗rightchoice)

	𝐔poststereoclick = poststereoclickbasis(Φpoststereoclick, trialdurations)
	𝐔premovement = premovementbasis(movementtimesteps, Φpremovement, trialdurations)
	𝐲 = vcat((trial.spiketrains[n] for trial in trials)...)
	Φdrift, 𝐔drift = drift_design_matrix(options, stereoclick_times_s, trialdurations, 𝐲)
	𝐔postspike = spikehistorybasis(Φpostspike, trialdurations, 𝐲)
	𝐗=hcat(𝐔drift, 𝐔postspike, 𝐔poststereoclick, 𝐔premovement, 𝐔postphotostimulus, 𝐕)
	indices𝐮 = Indices𝐮(size(𝐔drift,2), size(Φpostspike,2), size(Φpoststereoclick,2), size(Φpremovement,2), size(Φpostphotostimulus,2))
	glmθ = GLMθ(indices𝐮, size(𝐕,2), options)
	Model(options=options,
			Φdrift=Φdrift,
			Φpostspike=Φpostspike,
			Φpoststereoclick=Φpoststereoclick,
			Φpremovement=Φpremovement,
			𝐗=𝐗,
			𝐲=𝐲)
	end




	gaussianprior=GaussianPrior(options, trialsets)
	θnative = randomize_latent_parameters(options)
	θ₀native = FHMDDM.copy(θnative)
	Model(options=options,
		   gaussianprior=gaussianprior,
		   θnative=θnative,
		   θreal=native2real(options, θnative),
		   θ₀native=θ₀native,
		   trialsets=trialsets)
end
Model(options::Options) = Model(options, loadtrials(options))
