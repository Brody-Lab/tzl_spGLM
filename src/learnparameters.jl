"""
	maximizeposterior!(model)

Find the parameters of a Poisson GLM that gives the mode of the posterior distribution

MODIFIED ARGUMENT
-`model`: a struct containing the data, hyperparameters, and parameters. Only the parameters stored in the vector `𝐰` are modified
"""
function maximizeposterior!(model::Model; show_trace::Bool=true)
	x₀ = copy(model.𝐰)
	N = length(x₀)
	ℓ = fill(NaN,1)
	∇ℓ = fill(NaN, N)
	∇∇ℓ = fill(NaN, N, N)
	f(x) = neglogposterior!(model,ℓ,∇ℓ,∇∇ℓ,x)
	∇f!(∇,x) = ∇neglogposterior!(∇,model,ℓ,∇ℓ,∇∇ℓ,x)
	∇∇f!(∇∇,x) = ∇∇neglogposterior!(∇∇,model,ℓ,∇ℓ,∇∇ℓ,x)
    results = Optim.optimize(f, ∇f!, ∇∇f!, x₀, NewtonTrustRegion(), Optim.Options(show_trace=show_trace, iterations=model.options.opt_iterations_parameters))
	model.𝐰 .= Optim.minimizer(results)
	return nothing
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
function ∇∇neglogposterior!(∇∇::Matrix{<:Real}, model::Model, ℓ::Vector{<:Real}, ∇ℓ::Vector{<:Real}, ∇∇ℓ::Matrix{<:Real}, x::Vector{<:Real})
	∇∇logposterior!(model,ℓ,∇ℓ,∇∇ℓ,x)
	for i in eachindex(∇∇)
		∇∇[i] = -∇∇ℓ[i]
	end
	return nothing
end

"""
Gradient of the negative of the log-posterior
"""
function ∇neglogposterior!(∇::Vector{<:Real}, model::Model, ℓ::Vector{<:Real}, ∇ℓ::Vector{<:Real}, ∇∇ℓ::Matrix{<:Real}, x::Vector{<:Real})
	∇∇logposterior!(model,ℓ,∇ℓ,∇∇ℓ,x)
	for i in eachindex(∇)
		∇[i] = -∇ℓ[i]
	end
	return nothing
end

"""
Negative of the log-posterior
"""
function neglogposterior!(model::Model, ℓ::Vector{<:Real}, ∇ℓ::Vector{<:Real}, ∇∇ℓ::Matrix{<:Real}, x::Vector{<:Real})
	∇∇logposterior!(model,ℓ,∇ℓ,∇∇ℓ,x)
	-ℓ[1]
end

"""
Hessian of the log-posterior
"""
function ∇∇logposterior!(model::Model, ℓ::Vector{<:Real}, ∇ℓ::Vector{<:Real}, ∇∇ℓ::Matrix{<:Real}, x::Vector{<:Real})
	if (x != model.𝐰) || isnan(ℓ[1])
		for i in eachindex(x)
			model.𝐰[i] = x[i]
		end
		∇∇logposterior!(ℓ,∇ℓ,∇∇ℓ,model)
	end
	return nothing
end

"""
	∇∇logposterior!(ℓ, ∇ℓ, ∇∇ℓ, model)

Log-posterior and its gradient and hessian

MODIFIED ARGUMENT
-`ℓ`: log of the posterior distribution evaluated at `model.𝐰`
-`∇ℓ`: first-order derivatives of the expectation
-`∇∇ℓ`: second-order derivatives of the expectation

UNMODIFIED ARGUMENT
-`model`: a struct containing the data, parameters, and hyperparamters
"""
function ∇∇logposterior!(ℓ::Vector{<:Real}, ∇ℓ::Vector{<:Real}, ∇∇ℓ::Matrix{<:Real}, model::Model)
	∇∇loglikelihood!(ℓ, ∇ℓ, ∇∇ℓ, model)
    @unpack 𝐰 = model
	α = model.a[1]
	ℓ[1] -= α*(𝐰⋅𝐰)
	for i in eachindex(∇ℓ)
		∇ℓ[i] -= 2α*𝐰[i]
	    ∇∇ℓ[i,i] -= 2α
	end
	return nothing
end

"""
Log-likelihood and its gradient and hessian
"""
function ∇∇loglikelihood!(ℓ::Vector{<:Real}, ∇ℓ::Vector{<:Real}, ∇∇ℓ::Matrix{<:Real}, model::Model)
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
