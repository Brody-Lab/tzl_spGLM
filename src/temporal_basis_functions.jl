"""
	temporal_basis_functions(filtername, options)

Temporal basis functions used to fit a filter

ARGUMENT
-`filtername`: name of the filter
-`options`: settings of the model

RETURN
-`Φ`: temporal basis functions. Element Φ[τ,i] corresponds to the value of  i-th temporal basis function in the τ-th time bin in the kernel
"""
function temporal_basis_functions(basistype::String, options::Options)
	begin_s = getfield(options, Symbol("tbf_"*basistype*"_begin_s"))
	end_s = getfield(options, Symbol("tbf_"*basistype*"_end_s"))
	nvalues = ceil(Int, (end_s-begin_s) / options.dt)
	startingvalue = 1 + ((begin_s < 0.0) ? ceil(Int, -begin_s / options.dt) : 0)
	basisfunctions( nvalues,
					getfield(options, Symbol("tbf_"*basistype*"_D"));
					begins0=getfield(options, Symbol("tbf_"*basistype*"_begins0")),
					ends0=getfield(options, Symbol("tbf_"*basistype*"_ends0")),
                	η=getfield(options, Symbol("tbf_"*basistype*"_stretch")),
					startingvalue = startingvalue)
end

"""
    basisfunctions(nvalues, nfunctions)

Nonlinear basis functions of an one-dimensional integer input

ARGUMENT
-`nvalues`: number of values of the input tiled by the basis functions
-`N`: number of basis functions

OPTIONAL ARGUMENT
-`startingvalue`: an integer indicating the shift in the input values
-`η`: a positive scalar indicating the magnitude to which basis functions evaluated near the `startingvalue` is compressed and functions evaluated far from the `startingvalue` is stretched
-`begins0`: whether the value of each basis function at the first input value is 0
-`ends0`: whether the value of each basis function at the last input value is 0
-`orthogonal_to_ones:` whether the basis functions should be orthogonal to a constant vector

RETURN
-`Φ`: Matrix whose element Φ[i,j] corresponds to the value of the j-th temporal basis function at the i-th input
"""
function basisfunctions(nvalues::Integer, D::Integer; begins0::Bool=false, ends0::Bool=false, η::Real=0.0, startingvalue::Integer=1, orthogonal_to_ones::Bool=false, orthogonalize::Bool=true)
    if D == 1
        Φ = ones(nvalues,1)./nvalues
    else
        𝐱 = collect(1:nvalues) .- startingvalue
        if η > eps()
            𝐱 = asinh.(η.*𝐱)
        end
        Φ = raisedcosines(begins0, D, ends0, 𝐱)
        if orthogonal_to_ones
            Φ = orthogonalize_to_ones(Φ)
        end
        if orthogonalize
            Φ = orthonormalbasis(Φ)
        end
    end
    return Φ
end

"""
	temporal_basis_functions(begins0, Δt, duration_s, ends0, hz, scalefactor, stretch)

Value of each temporal basis at each time bin in a trialset

INPUT
-`begins0`: whether the basis begins at zero
-`Δt`: time bin, in seconds
-`duration_s`: duration in seconds
-`ends0`: whether the basis end at zero
-`hz`: number of temporal basis functions per second
-`scalefactor`: scaling
-`stretch`: nonlinear stretching of time

RETURN
-`Φ`: temporal basis functions. Element Φ[τ,i] corresponds to the value of  i-th temporal basis function in the τ-th time step in each trial
"""
function temporal_basis_functions(begins0::Bool, Δt::AbstractFloat, duration_s::Real, ends0::Bool, hz::Real, scalefactor::Real, stretch::Real; orthogonal_to_ones::Bool=false)
	nfunctions = ceil(hz*duration_s)
	if isnan(nfunctions) || (nfunctions < 1)
		return fill(1.0, 0, 0)
	else
		nfunctions = convert(Int, nfunctions)
		ntimesteps = ceil(Int, duration_s/Δt)
		temporal_basis_functions(begins0, ends0, nfunctions, ntimesteps, scalefactor, stretch; orthogonal_to_ones=orthogonal_to_ones)
	end
end

"""
	temporal_basis_functions(begins0, ends0, nfunctions, ntimesteps, scalefactor, stretch)

Value of each temporal basis at each time bin in a trialset

INPUT
-`begins0`: whether the basis begins at zero
-`Δt`: time bin, in seconds
-`ends0`: whether the basis end at zero
-`nfunctions`: number of temporal basis functions
-`ntimesteps`: number of time steps
-`scalefactor`: scaling
-`stretch`: nonlinear stretching of time

RETURN
-`Φ`: temporal basis functions. Element Φ[τ,i] corresponds to the value of  i-th temporal basis function in the τ-th time step in each trial
"""
function temporal_basis_functions(begins0::Bool, ends0::Bool, nfunctions::Integer, ntimesteps::Integer, scalefactor::Real, stretch::Real; orthogonal_to_ones::Bool=false)
	Φ = raisedcosines(begins0, ends0, nfunctions, ntimesteps, stretch)
	if orthogonal_to_ones
		Φ = orthogonalize_to_ones(Φ)
	end
	Φ = orthonormalbasis(Φ)
	Φ .*= scalefactor
	return Φ
end

"""
	temporal_basis_functions(Φ, 𝐓)

Value of each temporal basis at each time bin in a trialset

INPUT
-`Φ`: temporal basis functions. Element Φ[τ,i] corresponds to the value of  i-th temporal basis function in the τ-th time step in each trial
-`𝐓`: vector of the number of timesteps in each trial

RETURN
-`𝐕`: A matrix whose element 𝐕[t,i] indicates the value of the i-th temporal basis function in the t-th time bin in the trialset
"""
function temporal_basis_functions(Φ::Matrix{<:AbstractFloat}, 𝐓::Vector{<:Integer})
    Tmax = maximum(𝐓)
    D = size(Φ,2)
    𝐕 = zeros(sum(𝐓), D)
	if D > 0
	    k = 0
	    for T in 𝐓
	        for t = 1:T
	            k = k + 1;
	            𝐕[k,:] = Φ[t,:]
	        end
	    end
	end
	return 𝐕
end

"""
	premovementbasis(movementtimesteps, Φ, 𝐓)

Temporal basis functions for the premovement kernel

ARGUMENT
-`movementtimes_s`: time of movement relative to the stereoclick, in seconds
-`Φ`: temporal basis functions. Element Φ[τ,i] corresponds to the value of  i-th temporal basis function in the τ-th time bin in the kernel
-`𝐓`: number of timesteps

RETURN
-`𝐔`: A matrix whose element 𝐔[t,i] indicates the value of the i-th temporal basis function in the t-th time bin in the trialset
"""
function premovementbasis(movementtimesteps::Vector{<:Integer}, Φ::Matrix{<:AbstractFloat}, 𝐓::Vector{<:Integer})
	nbins, D = size(Φ)
	𝐔 = zeros(sum(𝐓), D)
	if D > 0
		τ = 0
		for i=1:length(𝐓)
			T = 𝐓[i]
			if movementtimesteps[i] < nbins
				j₀ = nbins - movementtimesteps[i] + 1
				for (t,j) in zip(1:T, j₀:nbins)
					𝐔[τ+t,:] = Φ[j,:]
				end
			else
				t₀ = movementtimesteps[i] - nbins + 1
				for (t,j) in zip(t₀:T, 1:nbins)
					𝐔[τ+t,:] = Φ[j,:]
				end
			end
			τ += T
		end
	end
	return 𝐔
end

"""
	spikehistorybasis(Φ, 𝐓, 𝐲)

Response of each temporal basis function parametrizing the postspike filter at each time step in the trialset

ARGUMENT
-`Φ`: Value of each temporal basis function parametrizing the postspike filter at each time step in the filter
-`𝐓`: number of time step in each trial for all trials in the trialset
-`𝐲`: number of spikes of one neuron in each time step, concatenated across trials

RETURN
-`𝐔`: a matrix whose element 𝐔[t,i] corresponds to the response of the i-th temporal basis function at the t-th time step in the trialset.
"""
function spikehistorybasis(Φ::Matrix{<:AbstractFloat}, 𝐓::Vector{<:Integer}, 𝐲::Vector{<:Integer})
	filterlength, D = size(Φ)
	𝐔 = zeros(sum(𝐓), D)
	if D > 0
		τ = 0
		for T in 𝐓
			indices = τ .+ (1:T)
			spiketimesteps = findall(𝐲[indices] .> 0)
			for tₛₚₖ in spiketimesteps
				y = 𝐲[tₛₚₖ+τ]
				indices𝐔 = τ .+ (tₛₚₖ+1:T)
				for (i,j) in zip(indices𝐔, 1:filterlength)
					for p = 1:D
						𝐔[i,p] += y*Φ[j,p]
					end
				end
			end
			τ = τ + T;
		end
	end
	return 𝐔
end

"""
	poststereoclickbasis(Φ, 𝐓)

Value of each temporal basis vector of the post-stereoclick filter at each time bin in a trialset

INPUT
-`Φ`: temporal basis functions. Element Φ[τ,i] corresponds to the value of  i-th temporal basis function in the τ-th time step in each trial
-`𝐓`: vector of the number of timesteps in each trial

RETURN
-`𝐔`: A matrix whose element 𝐔[t,i] indicates the value of the i-th temporal basis function in the t-th time bin in the trialset
"""
function poststereoclickbasis(Φ::Matrix{<:AbstractFloat}, 𝐓::Vector{<:Integer})
    ntimesteps, nfunctions = size(Φ)
    𝐔 = zeros(sum(𝐓), nfunctions)
    τ = 0
    for T in 𝐓
        for t = 1:min(T,ntimesteps)
			for i = 1:nfunctions
	            𝐔[τ+t,i] = Φ[t,i]
			end
        end
		τ+=T
    end
	return 𝐔
end

"""
	orthogonalize_to_ones(Φ)

Orthogonalize the columns of a matrix to a vector of ones

RETURN
-A matrix whose columns are orthogonal to any vector whose elements have the same value
"""
function orthogonalize_to_ones(Φ::Matrix{<:AbstractFloat})
	nrows = size(Φ,1)
	(I - fill(1/nrows,nrows,nrows))*Φ
end

"""
	orthonormalbasis(X)

Orthonormal basis for the real vector space spanned by the columns of `X`

OPTIONAL ARGUMENT
-`min_relative_singular_value`: dimensions whose singular value, relative to the maximum singular value across dimensions, is less than `min_relative_singular_value` are omitted

RETURN
-A unitary matrix whose columns span the same real vector space spanned by the columns of `X`
"""
function orthonormalbasis(Φ::Matrix{<:AbstractFloat}; min_relative_singular_value::AbstractFloat=1e-2)
	factorization = svd(Φ)
	relative_singular_values = factorization.S./maximum(factorization.S)
	indices = relative_singular_values .> min_relative_singular_value
	return factorization.U[:,indices]
end

"""
    raisedcosines(begins0, ends0, D, nbins, stretch)

Values of raised cosine temporal basis functions (tbf's)

ARGUMENT
-`begins0`: whether the first temporal basis function begins at the trough or at the peak
-`ends0`: whether the last temporal basis function begins at the trough or at the peak
-`D`: number of bases
-`nbins`: number of bins in the time window tiled by the bases
-`stretch`: an index of the stretching of the cosines

OPTIONAL ARGUMENT
-`period`: period of the raised cosine, in units of the inter-center distance

RETURN
-`Φ`: Matrix whose element Φ[i,j] corresponds to the value of the j-th temporal basis function at the i-th timestep from beginning of the trial
"""
function raisedcosines(begins0::Bool, ends0::Bool, D::Integer, nbins::Integer, stretch::Real; period::Real=4)
    if isnan(stretch) || stretch < eps()
        a = 1
        b = nbins
        t = collect(1:nbins)
    else
        λ = 1/stretch
        a = log(1+λ)
        b = log(nbins+λ)
        t = log.(collect(1:nbins) .+ λ)
    end
    if begins0
        if ends0
            Δcenter = (b-a) / (D+3)
        else
            Δcenter = (b-a) / (D+1)
        end
        centers = a .+ 2Δcenter .+ collect(0:max(1,D-1)).*Δcenter
    else
        if ends0
            Δcenter = (b-a) / (D+1)
        else
            Δcenter = (b-a) / (D-1)
        end
        centers = a .+ collect(0:max(1,D-1)).*Δcenter
    end
    ω = 2π/Δcenter/period
    Φ = raisedcosines(centers, ω, t)
    if !begins0 && !ends0 # allow temporal basis functions to parametrize a constant function
        lefttail = raisedcosines([centers[1]-Δcenter], ω, t)
        righttail = raisedcosines([centers[end]+Δcenter], ω, t)
        Φ[:,1] += lefttail
        Φ[:,end] += righttail
        indices = t .<= centers[1] + period/2*Δcenter
        deviations = 2.0 .- sum(Φ,dims=2) # introduced by time compression
        Φ[indices,1] .+= deviations[indices]
    end
    return Φ
end

"""
    raisedcosines(begins0, D, ends0, 𝐱)

Equally spaced raised cosine basis functions

ARGUMENT
-`begins0`: whether the first temporal basis function begins at the trough or at the peak
-`ends0`: whether the last temporal basis function begins at the trough or at the peak
-`D`: number of bases
-`𝐱`: input to the raised cosine function

RETURN
-`Φ`: Matrix whose element Φ[i,j] corresponds to the value of the j-th temporal basis function at the i-th input
"""
function raisedcosines(begins0::Bool, D::Integer, ends0::Bool, 𝐱::Vector{<:Real})
    Δx = 𝐱[end]-𝐱[1]
    if begins0
        if ends0
            Δcenter = Δx / (D+3)
        else
            Δcenter = Δx / (D+1)
        end
        centers = 𝐱[1] .+ 2Δcenter .+ collect(0:max(1,D-1)).*Δcenter
    else
        if ends0
            Δcenter = Δx / (D+1)
        else
            Δcenter = Δx / (D-1)
        end
        centers = 𝐱[1] .+ collect(0:max(1,D-1)).*Δcenter
    end
    ω = π/2/Δcenter
    Φ = raisedcosines(centers, ω, 𝐱)
    if !begins0 && !ends0 # allow temporal basis functions to parametrize a constant function
        lefttail = raisedcosines([centers[1]-Δcenter], ω, 𝐱)
        righttail = raisedcosines([centers[end]+Δcenter], ω, 𝐱)
        Φ[:,1] += lefttail
        Φ[:,end] += righttail
        indices = 𝐱 .<= centers[1] + 2/Δcenter
        deviations = 2.0 .- sum(Φ,dims=2) # introduced by time compression
        Φ[indices,1] .+= deviations[indices]
    end
    return Φ
end

"""
    raisedcosines(centers, ω, t)

Values of raised cosine temporal basis functions

ARGUMENT
-`centers`: Vector of the centers of the raised cosines
-`ω`: angular frequency
-`t`: values at which the temporal basis functions are evaluated

RETURN
-`Φ`: Matrix whose element Φ[i,j] corresponds to the value of the j-th temporal basis function at the i-th timestep from beginning of the trial
"""
function raisedcosines(centers::Vector{<:AbstractFloat}, ω::AbstractFloat, t::Vector{<:Real})
    T = t .- centers'
    (cos.(max.(-π, min.(π, ω.*T))) .+ 1)/2
end
