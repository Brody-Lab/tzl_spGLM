"""
	maximizeevidence!(model)

Learn both the parameters and hyperparameters by maximizing the approximate evidence

This optimization procedure alternately fixes the hyperparameters and learn the parameters by maximizing the posterior, and fixes the parameters and learn the hyperparameters by maximizing the evidence.

MODIFIED ARGUMENT
-`model`: structure containing the parameters, hyperparameters, and data. The parameters and hyperparameters are updated.

RETURN
-`eo`: a struct containing the optimization trace

OPTIONAL ARGUMENT
-`MAP_convergence_g_tol`: maximum norm of the gradient at the end of the MAP optimization for convergence to be considered. A lack of convergence to the mode of the posterior distribution indicates that the Laplace approximation of the posterior distribution is inaccurate and the log-evidence cannot be accurately approximated using the approach in this module.
-`verbose`: whether to display messages
"""
function maximizeevidence!(model::Model; verbose::Bool=true, MAP_convergence_g_tol::Real = 1e-4)
	eo = EvidenceOptimization(model)
	memory = MemoryForOptimization(model)
	for i = 1:model.options.opt_iterations_hyperparameters
	    Optim_results = maximizeposterior!(memory, model)
		eo.a[i] = model.a[1]
		eo.MAP_g_residual[i] = getfield(Optim_results, :g_residual)
		eo.𝐰[i] .= model.𝐰
		verbose && println("Evidence optimization iteration: ", i, ": precision (a) = ",  model.a[1])
		verbose && println("Evidence optimization iteration: ", i, ": norm of the residual gradient of the log-posterior (g_residual) = ",  eo.MAP_g_residual[i])
		if eo.MAP_g_residual[i] < MAP_convergence_g_tol
			verbose && println("Evidence optimization iteration: ", i, ": MAP optimization converged")
		else
			verbose && println("Evidence optimization iteration: ", i, ": MAP optimization did not converge, and therefore the optimization procedure is aborting.")
			break
		end
		∇∇loglikelihood!(memory, model)
		eo.hessian_loglikelihood[i] .= memory.∇∇ℓ
		eo.𝐸[i] = logevidence(memory, model)
		verbose && println("Evidence optimization iteration: ", i, ": approximate log-evidence (𝐸) = ", eo.𝐸[i])
		if i < model.options.opt_iterations_hyperparameters
			model.a[1] = maximizeevidence(memory, model)
		end
	end
	if maximum(eo.𝐸) > -Inf
		iteration = findmax(eo.𝐸)[2]
		model.𝐰 .= eo.𝐰[iteration]
		model.a[1] = eo.a[iteration]
	end
	return eo
end

"""
RETURN a struct containing the optimization trace
"""
function EvidenceOptimization(model::Model)
	N = length(model.𝐰)
	M = model.options.opt_iterations_hyperparameters
	𝐰 = collect(fill(NaN,N) for i = 1:M)
	EvidenceOptimization(a=fill(NaN,M),
						𝐸=fill(-Inf,M),
						hessian_loglikelihood = collect(fill(NaN,N,N) for i =1:M),
						MAP_g_residual=fill(NaN,M),
						𝐰 = collect(fill(NaN,N) for i = 1:M))
end

"""
	maximizeevidence(memory, model)

RETURN precision parameter that maximizes the approximate evidence

ARGUMENT
-`memory`: structure containing variables to be modified during computations
-`model`: structure containing the parameters, hyperparameters, and data.
"""
function maximizeevidence(memory::MemoryForOptimization, model::Model)
	f(x) = -logevidence(memory, model, x[1])
	results = optimize(f, exp.(model.a), LBFGS(); autodiff = :forward)
	exp(Optim.minimizer(results)[1])
end

"""
	logevidence(memory, model, x)

RETURNS the approximate log evidence

ARGUMENT
-`a`: precision parameter
"""
function logevidence(a::Real, memory::MemoryForOptimization, model::Model)
	𝐇 = memory.∇∇ℓ
	𝐰ₘₐₚ = model.𝐰
	a₀ = model.a[1]
	A = (a*I - 𝐇)
	B = ((a₀*I - 𝐇)*𝐰ₘₐₚ)
	if any(isnan(x) || isinf(x) for x in A) || any(isnan(x) || isinf(x) for x in B)
		-Inf
	else
		𝐰 = A \ B
		loglikelihood(model,𝐰) - 0.5a*(𝐰⋅𝐰) - 0.5logdet(I - 𝐇./a)
	end
end
logevidence(memory::MemoryForOptimization, model::Model) = logevidence(model.a[1], memory, model)
logevidence(memory::MemoryForOptimization, model::Model, x::Real) = logevidence(exp(x), memory, model)
