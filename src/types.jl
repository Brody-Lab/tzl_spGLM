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
	time_in_trial_begin_s::TF=-0.5
	"Time, in second, before which spikes are included on each trial, aligned to the reference event on that trial."
	time_in_trial_end_s::TF=2.0; @assert time_in_trial_end_s > time_in_trial_begin_s
	"""
	the hyperparameters below specify the parametrization of each set of basis functions
		- `begin_s`: time relative to the event at which the kernel begins
		- `end_s`: time relative to the event at which the kernel ends
		- `begins0`: whether the first value of the kernel begins at zero
		- `ends0`: whether the last value of the kernel begins at zero
		- `N`: whether the first value of the kernel begins at zero
		- `distortion`: nonlinearity in the parametrization of the kernel. Ignored unless larger than 0.
		- `distortion_s`: time when the distortion is maximized
	"""
	"click-aligned linear filter"
	bfs_click_begin_s::TF= 0.01
	bfs_click_begins0::TB=true
	bfs_click_D::TI=5
	bfs_click_distortion::TF=1.0
	bfs_click_distortion_s::TF=0.03
	bfs_click_end_s::TF=0.5; @assert bfs_click_end_s > bfs_click_begin_s
	bfs_click_ends0::TB=true
	"movement-aligned linear filter"
	bfs_movement_begin_s::TF = -1.0
	bfs_movement_end_s::TF = 0.5; @assert bfs_movement_end_s > bfs_movement_begin_s
	bfs_movement_begins0::TB=true
	bfs_movement_ends0::TB=false
	bfs_movement_D::TI=5
	bfs_movement_distortion::TF=0.2
	bfs_movement_distortion_s::TF=0.0
	"postspike filter"
	bfs_postspike_begin_s::TF=0.01
	bfs_postspike_end_s::TF=0.25; @assert bfs_postspike_end_s > bfs_postspike_begin_s
	bfs_postspike_begins0::TB=false
	bfs_postspike_ends0::TB=true
	bfs_postspike_D::TI=5
	bfs_postspike_distortion::TF=1.0
	bfs_postspike_distortion_s::TF=0.01
	"time in trial aligned to the reference event"
	bfs_time_in_trial_begin_s::TF=time_in_trial_begin_s
	bfs_time_in_trial_end_s::TF=time_in_trial_end_s
	bfs_time_in_trial_begins0::TB=false
	bfs_time_in_trial_ends0::TB=false
	bfs_time_in_trial_D::TI=5
	bfs_time_in_trial_distortion::TF=0.2
	bfs_time_in_trial_distortion_s::TF=0.0
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
	input_postspike::TB = true
	input_time_in_trial::TB = true; @assert input_time_in_trial
	"maximum number of iterations for learning the parameters"
	opt_iterations_parameters::TI = 20
	"maximum number of iterations for learning the hyperparameters "
	opt_iterations_hyperparameters::TI = 3
	"maximum number of iterations for learning the hyperparameters "
	opt_MAP_convergence_g_tol::TI = 3
	"absolute path of the folder where the model output, including the summary and predictions, are saved"
	outputpath::TS=""
	"event on each trial aligned to which spikes are counted"
	reference_event::TS="cpoke_in"
end

"""
    Trial

Sensory stimuli, behavior, and spike train on each trial
"""
@with_kw struct Trial{	TB<:Bool,
						TF<:AbstractFloat,
						TI<:Integer,
						TVI<:Vector{<:Integer},
						TARF<:AbstractRange{<:AbstractFloat}}
	"behavioral choice"
    choice::TB
	"source of each click: `0` indicates a left click, `1` a right click, and `2` a stereoclick"
	clicks_source::TVI
	"time step of each click"
	clicks_timestep::TVI
	"log of the ratio of the generative right and left click rate"
	γ::TF
	"time step of leaving the center port, relative to the time of the stereoclick, in seconds"
	movement_timestep::TI
	"time of the reference event"
	reference_time_s::TF
	"right edges of the time bins"
	timesteps_s::TARF
	"index of the trial in the trialset"
	trialindex::TI
	"spike trains"
	y::TVI
	"spike train before the trial used only as independent variable for estimating the post-spike filter and not as dependent variable"
	ypre::TVI
	"number of time steps"
	T::TI = length(timesteps_s)
end

"""
	BasisFunctionSet

A set of basis functions to parametrize a linear filter
"""
@with_kw struct BasisFunctionSet{MF<:Matrix{<:AbstractFloat}, S<:Symbol, TARF<:AbstractRange{<:AbstractFloat}}
	name::S
	"matrix containing the values of the basis functions. Columns are different basis functions, rows are different time steps"
	Φ::MF
	"time steps in second aligned to the reference event"
	timesteps_s::TARF
end

"""
	WeightIndices

Indices of the encoding weights of the inputs
"""
@with_kw struct WeightIndices{UI<:UnitRange{<:Integer}}
	click::UI
	leftclick::UI
	leftmovement::UI
	movement::UI
	postspike::UI
	rightclick::UI
	rightmovement::UI
	time_in_trial::UI
end

"""
    Model

Poisson generalized linear model
"""
@with_kw struct Model{TO<:Options,
					VF<:Vector{<:AbstractFloat},
					MF<:Matrix{<:AbstractFloat},
					VT<:Vector{<:Trial},
					VB<:Vector{<:BasisFunctionSet},
					WI<:WeightIndices,
					VI<:Vector{<:Integer}}
	"precision parameter of the Gaussian prior on weights"
	a::VF=rand(1)
    "fixed hyperparameters"
    options::TO
	"set of basis functions"
	basissets::VB
	"for identifying the column of the design matrix"
	weightindices::WI
	"trials"
	trials::VT
	"design matrix"
	𝐗::MF
	"Poisson observations"
	𝐲::VI
	"concatenated weights"
	𝐰::VF=rand(size(𝐗,2))
end

"""
	MemoryForOptimization
"""
@with_kw struct MemoryForOptimization{VR<:Vector{<:Real}, MR<:Matrix{<:Real}}
	"log-likelihood or log-posterior"
	ℓ::VR
	"gradient"
	∇ℓ::VR
	"hessian"
	∇∇ℓ::MR
end

"""
	EvidenceOptimization
"""
@with_kw struct EvidenceOptimization{VR<:Vector{<:Real}, VVR<:Vector{<:Vector{<:Real}}}
	"precision used for each iteration of MAP optimization"
	a::VR
	"approximate log-evidence evaluated at the end of each MAP optimization"
	𝐸::VR
	"Euclidean norm of gradient of the log posterior at the end of each MAP optimization"
	MAP_g_residual::VR
	"MAP parameters"
	𝐰::VVR
end
