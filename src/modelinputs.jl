"""
	inputs_each_timestep(bfs, inputname, trials)

RETURN columns of the design matrix corresponding to model input

ARGUMENT
-`bfs`: a set of basis functions
-`inputname`: name of the model input
-`trials`: vector of trials
"""
function inputs_each_timestep(bfs::BasisFunctionSet, inputname::Symbol, trials::Vector{<:Trial})
	if inputname == :movement
		movementinputs(bfs, NaN, trials)
	elseif inputname == :leftmovement
		movementinputs(bfs, false, trials)
	elseif inputname == :rightmovement
		movementinputs(bfs, true, trials)
	elseif inputname == :response
		responseinputs(bfs, NaN, trials)
	elseif inputname == :leftresponse
		responseinputs(bfs, false, trials)
	elseif inputname == :rightresponse
		responseinputs(bfs, true, trials)
	elseif inputname == :click
		clickinputs(bfs, NaN, trials)
	elseif inputname == :stereoclick
		clickinputs(bfs, 2, trials)
	elseif inputname == :leftclick
		clickinputs(bfs, 0, trials)
	elseif inputname == :rightclick
		clickinputs(bfs, 1, trials)
	elseif inputname == :postspike
		postspikeinputs(bfs, trials)
	elseif inputname == :time_in_trial
		vcat(collect(bfs.Φ[1:trial.T,:] for trial in trials)...)
	elseif inputname == :fixation
		fixationinputs(bfs, NaN, trials)
	else
		error("unrecognized input")
	end
end

"""
	clickinputs(bfs, laterality, trials)

RETURN columns of the design matrix containng time-varying inputs aligned to time of click onset

ARGUMENT
-`bfs`: a set of basis functions
-`laterality`: whether the movement is leftward (`0`), rightward (`1`), or non-laterality(`NaN`)
-`trials`: vector of trials
"""
function clickinputs(bfs::BasisFunctionSet, source::Real, trials::Vector{<:Trial})
	N, D = size(bfs.Φ)
	∑T = sum(collect(trial.T for trial in trials))
	𝐔 = zeros(∑T, D)
	dt = bfs.timesteps_s[2]-bfs.timesteps_s[1]
	latency = ceil(Int, bfs.timesteps_s[1]-dt)
	if D > 0
		τ = 0
		for trial in trials
			for (t, s) in zip(trial.clicks_timestep, trial.clicks_source)
				if isnan(source) || (source==s)
					t₀ = t + latency
					for (i,j) in zip(t₀:trial.T, 1:N)
						for d = 1:D
							𝐔[τ+i,d] += bfs.Φ[j,d]
						end
					end
				end
			end
			τ += trial.T
		end
	end
	return 𝐔
end

"""
	movementinputs(bfs, source, trials)

RETURN columns of the design matrix containng inputs related to time of movement onset

ARGUMENT
-`bfs`: a set of basis functions
-`laterality`: whether the movement is leftward (`0`), rightward (`1`), or non-laterality(`NaN`)
-`trials`: vector of trials
"""
function movementinputs(bfs::BasisFunctionSet, laterality::Real, trials::Vector{<:Trial})
	N, D = size(bfs.Φ)
	∑T = sum(collect(trial.T for trial in trials))
	𝐔 = zeros(∑T, D)
	if D > 0
		Na = sum(bfs.timesteps_s .<= 0.0)
		Nb = N - Na
		τ = 0
		for trial in trials
			if isnan(laterality) || (laterality==trial.choice)
				Ta = sum(trial.timesteps_s .< 0.0)
				Tb = trial.T - Ta
				if trial.movement_timestep >= Na
					t₀ = trial.movement_timestep-Na+1
					for (t,j) in zip(t₀:trial.T, 1:N)
						𝐔[τ+t,:] = bfs.Φ[j,:]
					end
				else
					j₀ = Na - trial.movement_timestep + 1
					for (t,j) in zip(1:trial.T, j₀:N)
						𝐔[τ+t,:] = bfs.Φ[j,:]
					end
				end
			end
			τ += trial.T
		end
	end
	return 𝐔
end

"""
	fixationinputs(bfs, source, trials)

RETURN columns of the design matrix containing inputs related to time of the onset fixation (i.e., 'cpoke_in')

ARGUMENT
-`bfs`: a set of basis functions
-`laterality`: whether the decision is leftward (`0`), rightward (`1`), or non-laterality(`NaN`)
-`trials`: vector of trials
"""
function fixationinputs(bfs::BasisFunctionSet, laterality::Real, trials::Vector{<:Trial})
	N, D = size(bfs.Φ)
	∑T = sum(collect(trial.T for trial in trials))
	𝐔 = zeros(∑T, D)
	dt = bfs.timesteps_s[2]-bfs.timesteps_s[1]
	bfs_timesteps = ceil.(Int, bfs.timesteps_s./dt)
	if D > 0
		τ = 0
		for trial in trials
			if isnan(laterality) || (laterality==trial.choice)
				j₀ = -trial.fixation_timestep
				for t = 1:trial.T
					j = j₀ + t
					if (j >= 1) || (j <= N)
						𝐔[τ+t,:] = bfs.Φ[j,:]
					end
				end
			end
			τ += trial.T
		end
	end
	return 𝐔
end

"""
	responseinputs(bfs, source, trials)

RETURN columns of the design matrix containng inputs related to time of response completion

ARGUMENT
-`bfs`: a set of basis functions
-`laterality`: whether the movement is leftward (`0`), rightward (`1`), or non-laterality(`NaN`)
-`trials`: vector of trials
"""
function responseinputs(bfs::BasisFunctionSet, laterality::Real, trials::Vector{<:Trial})
	N, D = size(bfs.Φ)
	∑T = sum(collect(trial.T for trial in trials))
	𝐔 = zeros(∑T, D)
	if D > 0
		Na = sum(bfs.timesteps_s .<= 0.0)
		Nb = N - Na
		τ = 0
		for trial in trials
			if isnan(laterality) || (laterality==trial.choice)
				Ta = sum(trial.timesteps_s .< 0.0)
				Tb = trial.T - Ta
				if trial.response_timestep >= Na
					t₀ = trial.response_timestep-Na+1
					for (t,j) in zip(t₀:trial.T, 1:N)
						𝐔[τ+t,:] = bfs.Φ[j,:]
					end
				else
					j₀ = Na - trial.response_timestep + 1
					for (t,j) in zip(1:trial.T, j₀:N)
						𝐔[τ+t,:] = bfs.Φ[j,:]
					end
				end
			end
			τ += trial.T
		end
	end
	return 𝐔
end

"""
	postspikeinputs(bfs, trials)

RETURN columns of the design matrix containng inputs related to spike history

ARGUMENT
-`bfs`: a set of basis functions
-`trials`: vector of trials
"""
function postspikeinputs(bfs::BasisFunctionSet, trials::Vector{<:Trial})
	N, D = size(bfs.Φ)
	∑T = sum(collect(trial.T for trial in trials))
	𝐔 = zeros(∑T, D)
	if D > 0
		dt = bfs.timesteps_s[2]-bfs.timesteps_s[1]
		latency = ceil(Int, bfs.timesteps_s[1]-dt)
		τ = 0
		for trial in trials
			for t in findall(trial.y .> 0)
				yₜ = trial.y[t]
				t₀ = t+latency
				for (i,j) in zip(t₀:trial.T, 1:N)
					for d = 1:D
						𝐔[τ+i,d] += yₜ*bfs.Φ[j,d]
					end
				end
			end
			Tpre = trial.Tpre
			for t in findall(trial.ypre .> 0)
				yₜ = trial.ypre[t]
				j₀ = Tpre-t+1
				for (i,j) in zip(1:Tpre, j₀:N)
					for d = 1:D
						𝐔[τ+i,d] += yₜ*bfs.Φ[j,d]
					end
				end
			end
			τ = τ + trial.T
		end
	end
	return 𝐔
end
