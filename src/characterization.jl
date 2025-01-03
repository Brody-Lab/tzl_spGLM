"""
	Characterization(model)

Compute quantities useful for characterizing the model

ARGUMENT
-a composite containing the data, parameters, and hyperparameter of a factoral hidden Markov drift-diffusion model

RETURN
-a composite containg quantities useful for understanding the model. See `types.jl` for details on each field of the composite.
"""

"""
	Characterization(testmodel, trainingmodel)
"""
function Characterization(testmodel::Model, trainingmodel::Model)
	𝐄𝐞 = SPGLM.externalinput(testmodel)
	inferredrate, autocorrelation = inferrate(𝐄𝐞, testmodel)
	Characterization(LL = loglikelihood_each_timestep(testmodel, trainingmodel),
					 externalinput = 𝐄𝐞,
					 inferredrate = inferredrate,
					 autocorrelation = autocorrelation,
					 observed_spiketrains = collect(trial.y for trial in testmodel.trials),
					 trialindices = collect(trial.trialindex for trial in testmodel.trials))
end
Characterization(model::Model) = Characterization(model,model)

"""
	loglikelihood_each_timestep(model, trials)

RETURN the log-likelihood on each time step given model parameters, relative to a homogeneous Poisson model

ARGUMENT
-`model`: struct containing parameters for evaluating the log-likelihood. The parameter of the time-homogeneous Poisson model is estimated from the data contained in this struct
-`trials`: vector of the struct `Trial`.

"""
function loglikelihood_each_timestep(testmodel::Model, trainingmodel::Model)
	𝐋 = testmodel.𝐗*testmodel.𝐰
	ℓs = collect((zeros(trial.T) for trial in testmodel.trials))
	λΔt_homo = mean(trainingmodel.𝐲)
	τ = 0
	for (ℓ,trial) in zip(ℓs,testmodel.trials)
		for t = 1:trial.T
			ℓ[t] = poissonloglikelihood(trainingmodel.options.dt, 𝐋[t+τ], trial.y[t]) - poissonloglikelihood(λΔt_homo, trial.y[t])
		end
		τ += trial.T
	end
	ℓs
end

"""
	externalinput(model)

RETURN a nested vector indicating the component of the pre-nonlinearity input external to the model. Element `𝐄𝐞sorted[i][t]` corresponds to time step t of trial i.

The component of the pre-activation input from the spikes within the trial is omitted. Length is the number of time steps summed across trials
"""
function externalinput(model::Model)
	indices = map(fieldnames(WeightIndices)) do fieldname
		if !(fieldname == :postspike)
			getfield(model.weightindices, fieldname)
		else
			1:0
		end
	end
	indices = vcat(indices...)
	𝐄𝐞 = view(model.𝐗, :, indices)*model.𝐰[indices]
	𝐡 = postspikefilter(model)
	if sum(𝐡)>0.0
		τ = 0
		for trial in model.trials
			Tpre = trial.Tpre
			for t in findall(trial.ypre .> 0)
				yₜ = trial.ypre[t]
				j₀ = Tpre-t+1
				for (i,j) in zip(1:Tpre, j₀:Tpre)
					𝐄𝐞[τ+i] += yₜ*𝐡[j]
				end
			end
			τ += trial.T
		end
	end
	𝐄𝐞sorted = collect(zeros(trial.T) for trial in model.trials)
	τ = 0
	for e in 𝐄𝐞sorted
		T = length(e)
		e .= 𝐄𝐞[τ .+ (1:T)]
		τ = τ + T
	end
	𝐄𝐞sorted
end

"""
	postspikefilter(model)

RETURN a vector representing the post-spike filter

The first element of the vector corresponds to the first time step after the spike.
"""
function postspikefilter(model::Model)
	index = findfirst((set)->set.name==:postspike, model.basissets)
	Tpre = length(model.basissets[index].timesteps_s)
	𝐡 = zeros(Tpre)
	if !isempty(model.weightindices.postspike)
		𝐡 += model.basissets[index].Φ*model.𝐰[model.weightindices.postspike]
	end
	𝐡
end

"""
	inferrate(model)

RETURN
-`𝚲`: a nested vector whose element `𝚲[i][t]` the firing rate on on time step `t` of trial `i`
-`𝐑`: a nested vector whose element `𝐑[i][t]` is the autocorrelation function for lag `t` on trial `i`

OPTIONAL ARGUMENT
-`nsamples`: number of samples to draw
"""
inferrate(model::Model; nsamples=model.options.sampling_N) = inferrate(externalinput(model), model;nsamples=nsamples)
function inferrate(𝐄𝐞::Vector{<:Vector{<:AbstractFloat}}, model::Model; nsamples=model.options.sampling_N)
	𝐡 = postspikefilter(model)
	𝚲 = collect((zeros(trial.T) for trial in model.trials))
	𝐑 = collect((zeros(trial.T-1) for trial in model.trials))
	τ = 0
	for (𝐄𝐞,𝛌,𝐫,trial) in zip(𝐄𝐞,𝚲,𝐑,model.trials)
		inferrate!(𝛌, 𝐫, model.options.dt, 𝐄𝐞, 𝐡, trial; nsamples=nsamples)
		τ += trial.T
	end
	𝚲, 𝐑
end

"""
	inferrate!(𝛌, Δt, 𝐄𝐞, 𝐡, trial)

Infer the rate on each time step of a trial

MODIFIED ARGUMENT
-`𝛌`: a vector of floats for in place computation

ARGUMENT
-`Δt`: duration of each time step, in seconds
-`𝐄𝐞`: input external to the model
-`𝐡`: post-spikefilter. First index indicates the first time step, plus latency, after a spike

OPTIONAL ARGUMENT
-`nsamples`: number of samples
"""
function inferrate!(𝛌::Vector{<:AbstractFloat}, 𝐫::Vector{<:AbstractFloat}, Δt::AbstractFloat, 𝐄𝐞::Vector{<:AbstractFloat}, 𝐡::Vector{<:AbstractFloat}, trial::Trial; nsamples::Integer=100)
	𝕪 = similar(trial.y)
	lags = 1:trial.T-1
	𝕣 = similar(𝐫)
	N = 0
	for s = 1:nsamples
		sample!(𝕪, Δt, 𝐄𝐞, 𝐡)
		for t = 1:trial.T
			𝛌[t] += 𝕪[t]/nsamples
		end
		if sum(𝕪)>0
			StatsBase.autocor!(𝕣,𝕪,lags)
			𝐫 .+= 𝕣
			N+=1
		end
	end
	𝐫 ./= N
end

"""
	sample!(a, c, 𝐄𝐞, 𝐡, mpGLM, 𝛚, 𝛕)

Generate a sample of spiking response on each time step of one trial

MODIFIED ARGUMENT
-`𝕪`: a vector of integers for in-place computation

ARGUMENT
-`Δt`: duration of each time step, in seconds
-`𝐄𝐞`: input external to the model
-`𝐡`: post-spikefilter. First index indicates the first time step, plus latency, after a spike
"""
function sample!(𝕪::Vector{<:Integer}, Δt::AbstractFloat, 𝐄𝐞::Vector{<:AbstractFloat}, 𝐡::Vector{<:AbstractFloat})
	max_spikes_per_step = floor(1000Δt)
	T = length(𝕪)
	Tpre = length(𝐡)
    for t = 1:T
		L = 𝐄𝐞[t]
		for lag = 1:min(Tpre, t-1)
			if 𝕪[t-lag] > 0
				L += 𝐡[lag]*𝕪[t-lag]
			end
		end
		λ = inverselink(L)
		𝕪[t] = min(rand(Poisson(λ*Δt)), max_spikes_per_step)
    end
	return nothing
end
