"""
	basis_function_sets()

RETURN a vector of symbols indicating the standard collection of basis sets

No ARGUMENT
"""
function basis_function_sets()
	[:click, :movement, :postspike, :response, :time_in_trial]
end

"""
	match_input_to_basis(inputname)

RETURN the name of a set of basis functions matching the name of an input

ARGUMENT: a Symbol indicating the name of an input
"""
function match_input_to_basis(inputname::Symbol)
	if (inputname == :movement) || (inputname == :leftmovement) || (inputname == :rightmovement)
		:movement
	elseif (inputname == :response) || (inputname == :leftresponse) ||(inputname == :rightresponse)
		:response
	elseif (inputname == :click) || (inputname == :leftclick) || (inputname == :rightclick) || (inputname == :stereoclick)
		:click
	elseif inputname == :postspike
		:postspike
	elseif inputname == :time_in_trial
		:time_in_trial
	else
		error("unrecognized name of input")
	end
end

"""
RETURN a set of basis functions matching the name of an input

ARGUMENT
-`basissets`: a set of basis functions
-`inputname`: a Symbol indicating the name of an input
"""
match_input_to_basis(basissets::Vector{<:BasisFunctionSet}, inputname::Symbol) = filter((basis)->match_input_to_basis(inputname)==basis.name, basissets)

"""
	BasisFunctionSet(setname, options)

RETURN a struct containing the values of a set of basis functions

ARGUMENT
-`setname`: a symbol identifying the basis function set
-`options`: a struct containing the fixed hyperparameters
"""
function BasisFunctionSet(setname::Symbol, options::Options)
	if setname==:time_in_trial
		begin_s = getfield(options, :time_in_trial_begin_s)
		end_s = getfield(options, :time_in_trial_end_s)
	else
		begin_s = getfield(options, Symbol("bfs_"*String(setname)*"_begin_s"))
		end_s = getfield(options, Symbol("bfs_"*String(setname)*"_end_s"))
	end
	D = getfield(options, Symbol("bfs_"*String(setname)*"_D"))
	begin_s = floor(begin_s/options.dt)*options.dt
	end_s = ceil(end_s/options.dt)*options.dt
	timesteps_s = (begin_s+options.dt):options.dt:end_s
	N = length(timesteps_s)
	distortion_s = getfield(options, Symbol("bfs_"*String(setname)*"_distortion_s"))
	ηindex = ceil(Int, (distortion_s-begin_s)/options.dt)
	Φ = basisfunctions(N, D;
					begins0=getfield(options, Symbol("bfs_"*String(setname)*"_begins0")),
					ends0=getfield(options, Symbol("bfs_"*String(setname)*"_ends0")),
					η=getfield(options, Symbol("bfs_"*String(setname)*"_distortion")),
					ηindex = ηindex)
	BasisFunctionSet(timesteps_s = timesteps_s,
					name=setname,
					Φ=Φ)
end

"""
    basisfunctions(nvalues, D)

Nonlinear basis functions of an one-dimensional integer input

ARGUMENT
-`nvalues`: number of values of the input tiled by the basis functions
-`N`: number of basis functions

OPTIONAL ARGUMENT
-`begins0`: whether the value of each basis function at the first input value is 0
-`ends0`: whether the value of each basis function at the last input value is 0
-`η` and `ηindex`: basis functions are compressed near `ηindex` and dilated far from `ηindex`. The magnitude of compression and dilation are quantified by `η`.
-`orthogonal_to_ones:` whether the basis functions should be orthogonal to a constant vector

RETURN
-`Φ`: Matrix whose element Φ[i,j] corresponds to the value of the j-th temporal basis function at the i-th input
"""
function basisfunctions(nvalues::Integer, D::Integer; begins0::Bool=false, ends0::Bool=false, η::Real=0.0, ηindex::Integer=1, orthogonal_to_ones::Bool=false, orthogonalize::Bool=true)
    if D == 1
        Φ = ones(nvalues,1)./nvalues
    else
        𝐱 = collect(1:nvalues) .- ηindex
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
    if !begins0 && !ends0 # allow temporal basis functions to parametrize a constant function of the form `f(t) = a*t` for scalar `a`
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
