"""
	maximizeevidence!(model)

Learn both the parameters and hyperparameters by maximizing the approximate evidence

This optimization procedure alternately fixes the hyperparameters and learn the parameters by maximizing the posterior, and fixes the parameters and learn the hyperparameters by maximizing the evidence.

MODIFIED ARGUMENT
-`model`: structure containing the parameters, hyperparameters, and data. The parameters and hyperparameters are updated.

OPTIONAL ARGUMENT
-`verbose`: whether to display messages
"""
function maximizeevidence!(model::Model; verbose::Bool=true, MAP_convergence_g_tol::Real = 1e-5)
	best𝐸 = -Inf
	best𝐰 = copy(model.𝐰)
	best𝑎 = model.a[1]
	memory = MemoryForOptimization(model)
	for i = 1:model.options.opt_iterations_hyperparameters
		verbose && println("Evidence optimization iteration: ", i, ": maximizing the log-posterior.")
	    Optim_results = maximizeposterior!(memory, model)
		MAP_converged = getfield(Optim_results, :g_residual) < MAP_convergence_g_tol
		if MAP_converged
			verbose && println("Evidence optimization iteration: ", i, ": MAP optimization converged")
		else
			verbose && println("Evidence optimization iteration: ", i, ": MAP optimization did not converge, and therefore the optimization procedure is aborting.")
			break
		end
		∇∇loglikelihood!(memory, model)
		𝐸 = logevidence(memory, model)
		if 𝐸 > best𝐸
			verbose && println("Evidence optimization iteration: ", i, ": the current log-evidence ( ", 𝐸, ") is greater than its previous value (", best𝐸, ").")
			best𝐸 = 𝐸
			best𝐰 .= model.𝐰
			best𝑎 = model.a[1]
		else
			verbose && println("Evidence optimization iteration: ", i, ": the current log-evidence ( ", 𝐸, ") is not greater than its previous value (", best𝐸, "), and therefore the optimization procedure is aborting.")
			break
		end
		if i==model.options.opt_iterations_hyperparameters
			verbose && println("Evidence optimization iteration ", i, ": the last iteration has been reached, and the optimization procedure is aborting.")
			break
		end
		anew = maximizeevidence(memory, model)
		verbose && println("Evidence optimization iteration ", i, ": precision (a) ", model.a[1], " → ", anew)
		model.a[1] = anew
	end
	if best𝐸 > -Inf
		model.𝐰 .= best𝐰
		model.a[1] = best𝑎[1]
	end
	return nothing
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
	𝐰 = (a*I - 𝐇) \ ((a₀*I - 𝐇)*𝐰ₘₐₚ)
	loglikelihood(model,𝐰) - 0.5a*(𝐰⋅𝐰) - 0.5logdet(I - 𝐇./a)
end
logevidence(memory::MemoryForOptimization, model::Model) = logevidence(model.a[1], memory, model)
logevidence(memory::MemoryForOptimization, model::Model, x::Real) = logevidence(exp(x), memory, model)
