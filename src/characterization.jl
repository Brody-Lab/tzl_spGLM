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
function Characterization(model::Model)
	Characterization(LL = loglikelihood_each_timestep(model),
					 inferredrate = inferrate(model))

end

"""
RETURN a nested vector whose element `ℓ[i][t]` the log-likelihood on on time step `t` of trial `i`
"""
function loglikelihood_each_timestep(model::Model)
	𝐋 = model.𝐗*model.𝐰
	ℓs = collect((zeros(trial.T) for trial in model.trials))
	τ = 0
	for (ℓ,trial) in zip(ℓs,model.trials)
		for t = 1:trial.T
			ℓ[t] = poissonloglikelihood(model.options.dt, 𝐋[t+τ], trial.y[t])
		end
		τ += trial.T
	end
	ℓs
end

"""
	externalinput(model)

RETURN a vector indicating the component of the pre-nonlinearity input external to the model

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
	𝐄𝐞
end

"""
	postspikefilter(model)

RETURN a vector representing the post-spike filter

The first element of the vector corresponds to the first time step after the spike.
"""
function postspikefilter(model::Model)
	Tpre = ceil(Int, (model.options.bfs_postspike_end_s-model.options.bfs_postspike_begin_s)/model.options.dt)
	𝐡 = zeros(Tpre)
	index = findall((set)->set.name==:postspike, model.basissets)
	if length(index)==1
		𝐡 += model.basissets[index[1]].Φ*model.𝐰[model.weightindices.postspike]
	end
	𝐡
end

"""
	inferrate(model)

RETURN a nested vector whose element `𝚲[i][t]` the log-likelihood on on time step `t` of trial `i`

OPTIONAL ARGUMENT
-`nsamples`: number of samples to draw
"""
function inferrate(model::Model; nsamples=model.options.sampling_N)
	𝐄𝐞 = externalinput(model)
	𝐡 = postspikefilter(model)
	𝚲 = collect((zeros(trial.T) for trial in model.trials))
	τ = 0
	for (𝛌,trial) in zip(𝚲,model.trials)
		inferrate!(𝛌, model.options.dt, 𝐄𝐞[τ+1:τ+trial.T], 𝐡, trial; nsamples=nsamples)
		τ += trial.T
	end
	𝚲
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
function inferrate!(𝛌::Vector{<:AbstractFloat}, Δt::AbstractFloat, 𝐄𝐞::Vector{<:AbstractFloat}, 𝐡::Vector{<:AbstractFloat}, trial::Trial; nsamples::Integer=100)
	𝕪 = similar(trial.y)
	for s = 1:nsamples
		sample!(𝕪, Δt, 𝐄𝐞, 𝐡)
		for t = 1:trial.T
			𝛌[t] += 𝕪[t]/nsamples
		end
	end
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

"""
	save(characterization, folderpath)

Save the characterizations of the model.

Each field of the composite `characterization` is saved within a separate file, with the same name as that of the field, within a folder whose absolute path is specified by `folderpath`.
"""
function save(characterization::Characterization, folderpath::String)
	if !isdir(folderpath)
		mkdir(folderpath)
		@assert isdir(folderpath)
	end
	for fieldname in fieldnames(Characterization)
		filepath = joinpath(folderpath, String(fieldname)*".mat")
		entry = getfield(characterization, fieldname)
		dict = Dict(String(fieldname)=>entry)
	    matwrite(filepath, dict)
	end
end
