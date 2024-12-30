"""
	perievent_time_histograms(𝚲,model)

RETURN a vector containing instances of type `PerieventTimeHistogram`

Each element of the vector is a struct containing a peri-event time histogram (PETH) aligned to a particular event in the trial, averaged across trials of a particular condition. The struct contains PETH computed from either observed spike trains or inferred spike rates

ARGUMENT
-`𝚲`: inferred spike rate, from the field `inferredrate` of the type `Characterization`
-`model`: a struct containing the data, hyperparameters, and parameters
"""
perievent_time_histograms(𝚲::Vector{<:Vector{<:AbstractFloat}}, model::Model) = perievent_time_histograms(model.basissets, 𝚲, model.options, model.trials)

function perievent_time_histograms(basissets::Vector{<:BasisFunctionSet}, 𝚲::Vector{<:Vector{<:AbstractFloat}}, options::Options, trials::Vector{<:Trial})
	reference_events = [options.reference_event; "movement"; "stereoclick"]
	for event in ["fixation"; "leftclick";  "rightclick"; "response"; "leftresponse"; "rightresponse"]
		getfield(options, Symbol("input_"*event)) && (reference_events = vcat(reference_events, event))
	end
	reference_events = unique(reference_events)
	peths = map(reference_events) do reference_event
				perievent_time_histograms(basissets,𝚲,options,reference_event,trials)
			end
	vcat(peths...)
end

"""
		perievent_time_histograms(basissets,𝚲,options,reference_event,trials)

RETURN a vector containing instances of type `PerieventTimeHistogram`, aligned to the String `reference_event`
"""
function perievent_time_histograms(basissets::Vector{<:BasisFunctionSet}, 𝚲::Vector{<:Vector{<:AbstractFloat}}, options::Options, reference_event::String, trials::Vector{<:Trial})
	Y = collect(trial.y for trial in trials)
	𝐲 = vcat(Y...)
	if reference_event==options.reference_event
		basisindex = first(findall((set)->set.name==:time_in_trial, basissets))
		if options.time_in_trial_end_s == 0.0
			reference_timesteps = collect([findfirst(trial.timesteps_s .== options.time_in_trial_begin_s)] for trial in trials)
		else
			reference_timesteps = collect([findfirst(trial.timesteps_s .> options.time_in_trial_begin_s)] for trial in trials)
		end
	elseif reference_event=="movement"
		basisindex = first(findall((set)->set.name==:movement, basissets))
		reference_timesteps = collect([trial.movement_timestep] for trial in trials)
	elseif reference_event=="fixation"
		basisindex = first(findall((set)->set.name==:fixation, basissets))
		reference_timesteps = collect([trial.fixation_timestep] for trial in trials)
	elseif reference_event=="response"
		basisindex = first(findall((set)->set.name==:response, basissets))
		reference_timesteps = collect([trial.response_timestep] for trial in trials)
	elseif reference_event=="postspike"
		basisindex = first(findall((set)->set.name==:postspike, basissets))
		ymax = maximum(𝐲)
		reference_timesteps = collect(vcat((repeat(findall(trial.y .== i),i) for i = 1:ymax)...) for trial in trials)
	elseif contains(reference_event, "click")
		basisindex = first(findall((set)->set.name==:click, basissets))
		if reference_event == "leftclick"
			source = 0.0
		elseif reference_event == "rightclick"
			source = 1.0
		elseif reference_event == "stereoclick"
			source = 2.0
		elseif reference_event == "click"
			source = NaN
		else
			error("unrecognized reference event")
		end
		reference_timesteps =
			map(trials) do trial
				if isnan(source)
					trial.clicks_timestep
				else
					trial.clicks_timestep[trial.clicks_source .== source]
				end
			end
	else
		error("unrecognized reference event")
	end
	timesteps_s = basissets[basisindex].timesteps_s
	timesteps = ceil.(Int, timesteps_s./options.dt)
	if reference_event==options.reference_event
		timesteps = collect(1:length(timesteps_s))
	else
		timesteps = ceil.(Int, timesteps_s./options.dt)
	end
	peths = map(collect(fieldnames(SPGLM.PETHSet))) do condition
				trialindices = collect(SPGLM.selecttrial(condition, trial) for trial in trials)
				if sum(trialindices)==0
					SPGLM.PerieventTimeHistogram(condition=String(condition),
											observed=fill(NaN, length(timesteps_s)),
											predicted=fill(NaN, length(timesteps_s)),
											reference_event=reference_event,
											timesteps_s=collect(timesteps_s))
				else
					observed = align_and_average(reference_timesteps[trialindices], timesteps, Y[trialindices])./options.dt
					predicted = align_and_average(reference_timesteps[trialindices], timesteps, 𝚲[trialindices])./options.dt
					SPGLM.PerieventTimeHistogram(condition=String(condition),
											observed=observed,
										   	predicted=predicted,
										   	reference_event=reference_event,
										   	timesteps_s=collect(timesteps_s))
				end
		end
end

"""
	align_and_average

Align time series and average

RETURN
-a vector of floats

ARGUMENT
-`reference_timesteps`: a nested vector of positive integers indicating the time step on which a reference event occurred in the source time series. Element `reference_timesteps[i][j]` indicates the time step when on trial `i`, the occurrence of the j-th event on that trial
-`timesteps`: a vector of positive integers indicating the time steps relative to the reference event of the PETH
-`X`: source time series. Element `X[i][t]` corresponds to time step `t` on trial `i`
"""
function align_and_average(reference_timesteps::Vector{<:Vector{<:Integer}}, timesteps::Vector{<:Integer}, X::Vector{<:Vector{<:Real}})
	N = length(timesteps)
	∑x = zeros(N)
	w = zeros(N)
	for (ts_ref,x) in zip(reference_timesteps,X)
		T = length(x)
		for tref in ts_ref
			for i = 1:N
				j = timesteps[i] + tref - 1
				if (j>=1) && (j<=T)
					∑x[i] += x[j]
					w[i] += 1
				end
			end
		end
	end
	∑x./w
end

"""
	selecttrial(condition, trial)

Does the trial fall under the given task condition?

ARGUMENT
-`condition`: a symbol naming the condition
-`trial`: a composite containing the data of a given trial

RETURN
-a Bool
"""
function selecttrial(condition::Symbol, trial::Trial)
	if condition==:unconditioned
		true
	elseif condition==:leftchoice
		!trial.choice
	elseif condition==:rightchoice
		trial.choice
	elseif condition==:leftevidence
		selecttrial(condition, trial.γ)
	elseif condition==:rightevidence
		selecttrial(condition, trial.γ)
	elseif condition==:leftchoice_strong_leftevidence
		(!trial.choice) && selecttrial(condition, trial.γ)
	elseif condition==:leftchoice_weak_leftevidence
		(!trial.choice) && selecttrial(condition, trial.γ)
	elseif condition==:rightchoice_strong_rightevidence
		(trial.choice) && selecttrial(condition, trial.γ)
	elseif condition==:rightchoice_weak_rightevidence
		(trial.choice) && selecttrial(condition, trial.γ)
	else
		error("unrecognized condition: "*String(condition))
	end
end

"""
	selecttrial(condition, γ)

Does the generative click rate fall under the given task condition?

ARGUMENT
-`condition`: a symbol naming the condition
-`γ`: log ratio of the generative right and left click rates

RETURN
-a Bool
"""
function selecttrial(condition::Symbol, γ::AbstractFloat)
	if condition==:leftevidence
		γ < 0
	elseif condition==:rightevidence
		γ > 0
	elseif condition==:leftchoice_strong_leftevidence
		γ < -2.25
	elseif condition==:leftchoice_weak_leftevidence
		(γ >= -2.25) && (γ < 0)
	elseif condition==:rightchoice_strong_rightevidence
		γ > 2.25
	elseif condition==:rightchoice_weak_rightevidence
		(γ <= 2.25) && (γ > 0)
	else
		error("unrecognized condition: "*String(condition))
	end
end

"""
	selecttrial(condition, Erightchoice, trial)

Does expected spike trains of a given trial fall under the given task condition?

ARGUMENT
-`condition`: a symbol naming the condition
-`Erightchoice`: expectation a random simulation of this trial results in a right choice
-`trial`: a composite containing the data of a given trial

RETURN
-a Bool
"""
function selecttrial(condition::Symbol, Erightchoice::AbstractFloat, trial::Trial)
	if condition==:unconditioned
		true
	elseif condition==:leftchoice
		Erightchoice < 1.0
	elseif condition==:rightchoice
		Erightchoice > 0.0
	elseif condition==:leftevidence
		selecttrial(condition, trial.γ)
	elseif condition==:rightevidence
		selecttrial(condition, trial.γ)
	elseif condition==:leftchoice_strong_leftevidence
		(Erightchoice < 1.0) && selecttrial(condition, trial.γ)
	elseif condition==:leftchoice_weak_leftevidence
		(Erightchoice < 1.0) && selecttrial(condition, trial.γ)
	elseif condition==:rightchoice_strong_rightevidence
		(Erightchoice > 0.0) && selecttrial(condition, trial.γ)
	elseif condition==:rightchoice_weak_rightevidence
		(Erightchoice > 0.0) && selecttrial(condition, trial.γ)
	else
		error("unrecognized condition: "*String(condition))
	end
end
