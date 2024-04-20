"""
	maximizeposterior!(memory, model)

MODIFIED ARGUMENT
-`memory`: a struct pointing to pre-allocated memory for model optimization
-`model`: a struct containing the data, hyperparameters, and parameters. Only the parameters stored in the vector `𝐰` are modified

RETURN
-a struct containing the results of the optimization (an instance of the composite type `Optim.MultivariateOptimizationResults`)
"""
function maximizeposterior!(memory::MemoryForOptimization, model::Model)
	memory.ℓ[1] = NaN
	f(x) = neglogposterior!(memory,model,x)
	∇f!(∇,x) = ∇neglogposterior!(∇,memory,model,x)
	∇∇f!(∇∇,x) = ∇∇neglogposterior!(∇∇,memory,model,x)
    results = Optim.optimize(f, ∇f!, ∇∇f!, copy(model.𝐰), NewtonTrustRegion(), Optim.Options(show_trace=true, iterations=model.options.opt_iterations_parameters))
	model.𝐰 .= Optim.minimizer(results)
	return results
end
maximizeposterior!(model::Model) = maximizeposterior!(MemoryForOptimization(model),model)

"""
	MemoryForOptimization(model)

RETURN a struct pointing to pre-allocated memory for model optimization
"""
function MemoryForOptimization(model::Model)
	N = length(model.𝐰)
	MemoryForOptimization(ℓ = fill(NaN,1),
						  ∇ℓ = fill(NaN, N),
						  ∇∇ℓ = fill(NaN, N, N))
end

"""
	∇∇neglogposterior!(∇∇,model,ℓ,∇ℓ,∇∇ℓ,x)

Hessian of the negative of the log-posterior

MODIFIED ARGUMENT
-`∇∇`: hessian matrix
-`model`: the field `𝐰` is modified
-`ℓ`: log-posterior
-`∇ℓ`: gradient of the log-posterior
-`∇∇ℓ`: hessian of the log-posterior

UNMODIFIED ARGUMENT
-`x`: parameter values
"""
function ∇∇neglogposterior!(∇∇::Matrix{<:Real}, memory::MemoryForOptimization, model::Model, x::Vector{<:Real})
	∇∇logposterior!(memory,model,x)
	for i in eachindex(∇∇)
		∇∇[i] = -memory.∇∇ℓ[i]
	end
	return nothing
end

"""
Gradient of the negative of the log-posterior
"""
function ∇neglogposterior!(∇::Vector{<:Real}, memory::MemoryForOptimization, model::Model, x::Vector{<:Real})
	∇∇logposterior!(memory,model,x)
	for i in eachindex(∇)
		∇[i] = -memory.∇ℓ[i]
	end
	return nothing
end

"""
Negative of the log-posterior
"""
function neglogposterior!(memory::MemoryForOptimization, model::Model, x::Vector{<:Real})
	∇∇logposterior!(memory,model,x)
	-memory.ℓ[1]
end

"""
Hessian of the log-posterior
"""
function ∇∇logposterior!(memory::MemoryForOptimization, model::Model, x::Vector{<:Real})
	if (x != model.𝐰) || isnan(memory.ℓ[1])
		for i in eachindex(x)
			model.𝐰[i] = x[i]
		end
		∇∇logposterior!(memory,model)
	end
	return nothing
end

"""
	∇∇logposterior!(memory, model)

Log-posterior and its gradient and hessian

MODIFIED ARGUMENT
-`memory`: a struct pointing to pre-allocated memory for model optimization

UNMODIFIED ARGUMENT
-`model`: a struct containing the data, parameters, and hyperparamters
"""
function ∇∇logposterior!(memory::MemoryForOptimization, model::Model)
	∇∇loglikelihood!(memory, model)
	@unpack ℓ, ∇ℓ, ∇∇ℓ = memory
    @unpack 𝐰 = model
	α = model.a[1]
	ℓ[1] -= 0.5α*(𝐰⋅𝐰)
	for i in eachindex(∇ℓ)
		∇ℓ[i] -= α*𝐰[i]
	    ∇∇ℓ[i,i] -= α
	end
	return nothing
end

"""
Log-likelihood and its gradient and hessian
"""
function ∇∇loglikelihood!(memory::MemoryForOptimization, model::Model)
	@unpack ℓ, ∇ℓ, ∇∇ℓ = memory
    @unpack 𝐗, 𝐰, 𝐲 = model
	Δt = model.options.dt
	𝐋 = 𝐗*𝐰
    ℓ[1] = 0.0
    d²𝐥_d𝐋², d𝐥_d𝐋 = similar(𝐋), similar(𝐋)
    for t in eachindex(𝐋)
        d²𝐥_d𝐋²[t], d𝐥_d𝐋[t], ℓₜ = differentiate_twice_loglikelihood_wrt_linearpredictor(Δt, 𝐋[t], 𝐲[t])
	    ℓ[1] += ℓₜ
    end
    ∇ℓ .= 𝐗'*d𝐥_d𝐋
    ∇∇ℓ .= 𝐗'*(d²𝐥_d𝐋².*𝐗)
	return nothing
end

"""
Function used for testing with machine differentiation
"""
logposterior(model::Model, 𝐰::Vector{<:Real}) = loglikelihood(model, 𝐰) - model.a[1]*(𝐰⋅𝐰)

"""
Function used for testing with machine differentiation
"""
function loglikelihood(model::Model, 𝐰::Vector{<:Real})
	@unpack 𝐗, 𝐲 = model
	Δt = model.options.dt
	𝐋 = 𝐗*𝐰
	ℓ = 0
	for (L,y) in zip(𝐋,𝐲)
		ℓ += poissonloglikelihood(Δt, L, y)
	end
	return ℓ
end
