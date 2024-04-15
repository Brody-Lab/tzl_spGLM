"""
    Options

Model settings
"""
@with_kw struct Options{TB<:Bool, TF<:AbstractFloat, TI<:Integer, TS<:String}
	"absolute path of the file containing the data"
    datapath::TS=""
	"duration of each timestep in seconds"
    dt::TF=0.01
	"Time, in second, after which spikes are included on each trial, aligned to the reference event on that trial."
	reference_begin_s::TF=-0.5
	"Time, in second, before which spikes are included on each trial, aligned to the reference event on that trial."
	reference_end_s::TF=2.0
	"""
	the hyperparameters below specify the parametrization of the linear kernels
		- `begin_s`: time relative to the event at which the kernel begins
		- `end_s`: time relative to the event at which the kernel ends
		- `begins0`: whether the first value of the kernel begins at zero
		- `ends0`: whether the last value of the kernel begins at zero
		- `N`: whether the first value of the kernel begins at zero
		- `stretch`: nonlinearity in the parametrization of the kernel
	"""
	"click-aligned linear filter"
	tbf_click_begin_s::TF= -0.01
	tbf_click_end_s::TF=0.49
	tbf_click_begins0::TB=true
	tbf_click_ends0::TB=true
	tbf_click_D::TI=5
	tbf_click_stretch::TF=1.0
	"movement-aligned linear filter"
	tbf_movement_begin_s::TF = -1.0
	tbf_movement_end_s::TF = 0.5
	tbf_movement_begins0::TB=true
	tbf_movement_ends0::TB=false
	tbf_movement_D::TI=5
	tbf_movement_stretch::TF=0.2
	"postspike filter"
	tbf_postspike_begin_s::TF=0.0
	tbf_postspike_end_s::TF=0.25
	tbf_postspike_begins0::TB=false
	tbf_postspike_ends0::TB=true
	tbf_postspike_D::TI=5
	tbf_postspike_stretch::TF=1.0
	"time in trial aligned to the reference event"
	tbf_reference_begin_s::TF=reference_begin_s
	tbf_reference_end_s::TF=reference_end_s
	tbf_reference_begins0::TB=false
	tbf_postspike_ends0::TB=false
	tbf_reference_D::TI=5
	tbf_reference_stretch::TF=1.0
	"""
	hyperparameters governing the inputs
	"""
	input_click::TB = false
	input_leftclick::TB = true
	input_rightclick::TB = true
	input_stereoclick::TB = false
	input_movement::TB = false
	input_leftmovement::TB = true
	input_rightmovement::TB = true
	"latency of click-related input"
	latency_click_s::TF = 0.01
	"latency of the post-spike input"
	latency_postspike_s::TF = 0.01
	"absolute path of the folder where the model output, including the summary and predictions, are saved"
	outputpath::TS=""
	"event on each trial aligned to which spikes are counted"
	reference_event::TS="cpoke_in"
end

"""
    Trial

Information on the sensory stimulus and behavior each trial

Spike trains are not included. In sampled data, the generatives values of the latent variables are stored.
"""
@with_kw struct Trial{	TB<:Bool,
						TF<:AbstractFloat,
						TVF<:Vector{<:AbstractFloat},
						TI<:Integer,
						TVI<:Vector{<:Integer},
						TARF<:AbstractRange{<:AbstractFloat}}
	"edges of the time bins"
	binedges_s::TARF
	"behavioral choice"
    choice::TB
	"time of the onset of nose fixation"
	fixation_time_s::TF
	"log of the ratio of the generative right and left click rate"
	γ::TF
	"times of left clicks"
	Lclick_times_s::TVF
	"time of leaving the center port, relative to the time of the stereoclick, in seconds"
	movement_time_s::TF
	"times of right clicks"
	Rclick_times_s::TVF
	"time of the reference event"
	reference_time_s::TF
	"time of the stereoclick, in seconds, in the sessions"
	stereoclick_time_s::TF
	"index of the trial in the trialset"
	trialindex::TI
	"spike trains"
	y::TVI
end

"""
	WeightIndices

Indices of the encoding weights of the inputs
"""
@with_kw struct WeightIndices{UI<:UnitRange{<:Integer}}
	click::UI
	fixation::UI
	leftclick::UI
	leftmovement::UI
	movement::UI
	postspike::UI
	reference::UI
	rightclick::UI
	rightmovement::UI
	stereoclick::UI
end

"""
	GLMKernel
"""
@with_kw struct GLMKernel{TB<:Bool, TF<:AbstractFloat, TI<:Integer, TS<:Symbol, MF<:Matrix{<:AbstractFloat}}
	"name"
	name::TS
	"columns correspond to temporal basis functions, rows to the value of a function at each time step"
	Φ::MF
	"start time relative to the reference event"
	begin_s::TF
	"end time relative to the reference event"
	end_s::TF
	"is the value of the filter equal to 0 at the beginning"
	begin_0::TB
	"is the value of the filter equal to 0 at the end"
	end_0::TB
	"number of temporal basis functions"
	N::TI
	"stretch in the nonlinearity"
	η::TF
end

"""
	GaussianPrior

Information on the zero-meaned Gaussian prior distribution on the values of the parameters in real space
"""
@with_kw struct GaussianPrior{	VI<:Vector{<:Integer},
								VF<:Vector{<:AbstractFloat},
								VR<:Vector{<:Real},
								VS<:Vector{<:String},
								MR<:Matrix{<:Real},
								VVI<:Vector{<:Vector{<:Integer}},
								VMF<:Vector{<:Matrix{<:AbstractFloat}}}
	"L2 penalty matrices"
	𝐀::VMF
	"L2 penalty coefficients"
	𝛂::VR
	"minimum values of the L2 penalty coefficients"
	𝛂min::VF
	"maximum values of the L2 penalty coefficients"
	𝛂max::VF
	"Indices of the parameters related to each L2 penalty coefficient: element `index𝐀[i][j]` corresponds to the i-th group of parameters and the j-th parameter in that group"
	index𝐀::VVI
	"the precision matrix, i.e., inverse of the covariance matrix, of the gaussian prior on the model parameters"
	𝚲::MR
	"indices of the dimensions with finite variance"
	index𝚽::VI = sort(union(index𝐀...))
	"square submatrix of the precision matrix after deleting the columns and rows corresponding to the dimensions with infinite variance"
	𝚽::MR= 𝚲[index𝚽,index𝚽]
	"indices of 𝐀 within `index𝚽`"
	index𝐀_in_index𝚽::VVI = map(indexA->map(indexAᵢⱼ->findfirst(index𝚽.==indexAᵢⱼ), indexA), index𝐀)
	"the name of each L2 penalty"
	penaltynames::VS
end

"""
    Model

Poisson generalized linear model
"""
@with_kw struct Model{TO<:Options,
					VF<:Vector{<:AbstractFloat},
					MF<:Matrix{<:AbstractFloat},
					VK<:Vector{<:GLMKernel},
					VT<:Vector{<:Trial},
					GP<:GaussianPrior,
					WI<:WeightIndices,
					VI<:Vector{<:UInt8}}
	"gaussian prior"
	gaussianprior::GP
    "fixed hyperparameters"
    options::TO
	"linear filters"
	glmkernels::VK
	""
	weightindices::WI
	"trials"
	trials::VT
	"concatenated weights"
	𝐰::VF
	"design matrix"
	𝐗::MF
	"Poisson observations"
	𝐲::VI
end
