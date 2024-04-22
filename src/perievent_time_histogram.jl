"""
	perievent_time_histograms(𝚲,model)

RETURN a vector containing instances of type `PerieventTimeHistogram`

Each element of the vector is a struct containing a peri-event time histogram (PETH) aligned to a particular event in the trial, averaged across trials of a particular condition. The struct contains PETH computed from either observed spike trains or inferred spike rates

ARGUMENT
-`𝚲`: inferred spike rate, from the field `inferredrate` of the type `Characterization`
-`model`: a struct containing the data, hyperparameters, and parameters
"""
function perievent_time_histograms(𝚲::Vector{<:Vector{<:AbstractFloat}}, model::Model)
	reference_events = [model.options.reference_event; "movement"; "postspike"; "leftclick"; "rightclick"; "stereoclick"; "click"]
	peths = map(reference_events) do reference_event
				perievent_time_histograms(𝚲,model,reference_event)
			end
	vcat(peths...)
end

"""
		perievent_time_histograms(𝚲,model,reference_event)

RETURN a vector containing instances of type `PerieventTimeHistogram`, aligned to the String `reference_event`
"""
function perievent_time_histograms(𝚲::Vector{<:Vector{<:AbstractFloat}}, model::Model, reference_event::String)
	if reference_event==model.options.reference_event
		basisindex = first(findall((set)->set.name==:time_in_trial, model.basissets))
		reference_timesteps = collect([findfirst(trial.timesteps_s .> 0)] for trial in model.trials)
	elseif reference_event=="movement"
		basisindex = first(findall((set)->set.name==:movement, model.basissets))
		reference_timesteps = collect([trial.movement_timestep] for trial in model.trials)
	elseif reference_event=="postspike"
		basisindex = first(findall((set)->set.name==:postspike, model.basissets))
		ymax = maximum(model.𝐲)
		reference_timesteps = collect(vcat((repeat(findall(trial.y .== i),i) for i = 1:ymax)...) for trial in model.trials)
	elseif contains(reference_event, "click")
		basisindex = first(findall((set)->set.name==:click, model.basissets))
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
			map(model.trials) do trial
				if isnan(source)
					trial.clicks_timestep
				else
					trial.clicks_timestep[trial.clicks_source .== source]
				end
			end
	else
		error("unrecognized reference event")
	end
	timesteps_s = model.basissets[basisindex].timesteps_s
	Na = findfirst(timesteps_s.>0)
	Nb = length(timesteps_s)-Na
	Y = collect(trial.y for trial in model.trials)
	peths = map(collect(fieldnames(SPGLM.PETHSet))) do condition
		trialindices = collect(SPGLM.selecttrial(condition, trial) for trial in model.trials)
		observed = SPGLM.align_and_average(reference_timesteps[trialindices], Na, Nb, Y[trialindices])
		predicted = SPGLM.align_and_average(reference_timesteps[trialindices], Na, Nb, 𝚲[trialindices])
		SPGLM.PerieventTimeHistogram(condition=String(condition),
								observed=observed,
							   	predicted=predicted,
							   	reference_event=reference_event,
							   	timesteps_s=collect(timesteps_s))
		end
end

"""
	align_and_average

Align time series and average

ARGUMENT
-`reference_timesteps`: a nested vector of positive integers indicating the time step on which a reference event occurred in the source time series. Element `reference_timesteps[i][j]` indicates the time step when on trial `i`, the occurrence of the j-th event on that trial
-`Na`: number of time steps before and including the reference event in the target time series
-`Nb`: number of time steps after the reference event in the target time series
-`X`: source time series. Element `X[i][t]` corresponds to time step `t` on trial `i`
"""
function align_and_average(reference_timesteps::Vector{<:Vector{<:Integer}}, Na::Integer, Nb::Integer, X::Vector{<:Vector{<:Real}})
	N = Na + Nb
	∑x = zeros(N)
	w = zeros(N)
	for (ts_ref,x) in zip(reference_timesteps,X)
		for tref in ts_ref
			tb = min(length(x)-tref, Nb)
			ib = tref + tb
			jb = Na + tb
			ta = min(tref-1, Na-1)
			ia = tref - ta
			ja = Na - ta
			if (ia>0) & (ib>0) & (ja>0) & (jb>0)
				for (i,j) in zip(ia:ib, ja:jb)
					∑x[j] += x[i]
					w[j] +=1
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
