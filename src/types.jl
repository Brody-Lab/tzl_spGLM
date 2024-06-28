"""
    Options

Model settings
"""
@with_kw struct Options{TB<:Bool, TF<:AbstractFloat, TI<:Integer, TS<:String}
	"""
	hyperparameters for estimating the component of the baseline from changes across trials
	"""
	baseline_L2max::TF=1e10
	baseline_L2min::TF=1e1
	baseline_L2n::TI=10
	"absolute path to a `.MAT` file containing spike counts from a time period before the trial to be used to estimate the across-trial contribution to the baseline. This file should contain a matrix of floats of size #trials-by-#neurons. All neurons should be simultaneously recorded, and includes the spike counts from the neuron to which the GLM is being fitted as well."
	baseline_pretrial_spikecounts_path::TS=""
	"whether to use time in session to estimate the across-trial contribution to the baseline"
	baseline_time_in_session::TB=false
	"if using time in session for quantifying the trial-varying component of the baseline, this parameter specifies the number of basis functions to use"
	baseline_time_in_session_D::TI=8
	"""
	the hyperparameters below specify the parametrization of each set of basis function ("bfs"). Each hyparameter is associated with an event in the Poisson Clicks task, such as a click, movement (i.e., "cpoke_out"), or response (i.e., "spoke"), or associated with the neuron's own spiking (i.e., spike history). The hyperparameters attributed to "time_in_trial" are related to the whichever task event set as the reference event to which spike trains are aligned on each trial. This is typically the stereoclick.
		- `begin_s`: time relative to the event at which the kernel begins
		- `begins0`: whether the filter must begin with the value 0
		- `D`: number of basis functions
		- `distortion`: the magnitude to which basis functions are compressed close in time to the event and dilated far in time from the event.
		- `distortion_s`: time of maximal distortion. If this is set to 0, then maximal distortion occurs at the time of the event.
		- `end_s`: time relative to the event at which the kernel ends
		- `ends0`: whether the filter must end with the value 0
	"""
	"linear filter aligned to each click (left, right, and stereo)"
	bfs_click_begin_s::TF= 0.01
	bfs_click_begins0::TB=true
	bfs_click_D::TI=3
	bfs_click_distortion::TF=0.1
	bfs_click_distortion_s::TF=0.05
	bfs_click_end_s::TF=0.3; @assert bfs_click_end_s > bfs_click_begin_s
	bfs_click_ends0::TB=false
	"linear filter aligned to the movement (i.e., 'cpoke_out')"
	bfs_movement_begin_s::TF = -1.0
	bfs_movement_end_s::TF = -0.01; @assert bfs_movement_end_s > bfs_movement_begin_s
	bfs_movement_begins0::TB=true
	bfs_movement_ends0::TB=false
	bfs_movement_D::TI=3
	bfs_movement_distortion::TF=0.1
	bfs_movement_distortion_s::TF=-0.01
	"linear filter aligned to the response (i.e., 'spoke')"
	bfs_response_begin_s::TF = -0.5
	bfs_response_end_s::TF = 0.01; @assert bfs_response_end_s > bfs_response_begin_s
	bfs_response_begins0::TB=true
	bfs_response_ends0::TB=false
	bfs_response_D::TI=4
	bfs_response_distortion::TF=0.0
	bfs_response_distortion_s::TF=0.0
	"postspike filter"
	bfs_postspike_begin_s::TF=0.01
	bfs_postspike_end_s::TF=0.25; @assert bfs_postspike_end_s > bfs_postspike_begin_s
	bfs_postspike_begins0::TB=false
	bfs_postspike_ends0::TB=false
	bfs_postspike_D::TI=3
	bfs_postspike_distortion::TF=1.0
	bfs_postspike_distortion_s::TF=0.01
	"time in trial aligned to the reference event"
	bfs_time_in_trial_begins0::TB=false
	bfs_time_in_trial_ends0::TB=false
	bfs_time_in_trial_D::TI=3
	bfs_time_in_trial_distortion::TF=0.1
	bfs_time_in_trial_distortion_s::TF=0.0
	"absolute path of the file containing the data"
	datapath::TS=""
	"duration of each timestep in seconds"
	dt::TF=0.01
	"""
	hyperparameters specifying whether include an inputs. Be careful to avoid redundancies, such as setting both "input_click=true" and "input_leftclick=true," which contributes to trade-offs in the parameter estimates
	"""
	"all clicks: left, right, and stereo"
	input_click::TB = false
	input_leftclick::TB = true
	input_rightclick::TB = true
	input_stereoclick::TB = true
	input_movement::TB = false
	"trial ending in a left choice only"
	input_leftmovement::TB = true
	input_rightmovement::TB = true
	"whether `'`pose summary` acts as an input. A file path must be provided as `pose_filepath`"
	input_pose::TB = false
	input_postspike::TB = true
	input_response::TB = false
	input_leftresponse::TB = false
	input_rightresponse::TB = false
	input_time_in_trial::TB = true; @assert input_time_in_trial
	"maximum number of iterations for learning the parameters"
	opt_iterations_parameters::TI = 20
	"maximum number of iterations for learning a single hyperparameter that serves as the precision of the Gaussian prior on the encoding weights"
	opt_iterations_hyperparameters::TI = 3
	"optimization method"
	opt_method::TS="evidenceoptimization"
	"absolute path to a file containing the pose measurements. The file should be a MAT file containing an variable `posesignals` that is a T-by-D matrix, where T is the number of samples and D the dimensionality of the signal. The file should also contain a variable `firstframetime_s` that indicates the time of the first frame and a variable `samplingrate_hz` indicating the sampling rate in hertz"
	pose_filepath::TS=""
	"default value of the precision hyperparameter"
	precision::TF=1e-2
	"number of samples used to compute the expected emissions and peri-event time histograms"
	sampling_N::TI = 100
	"absolute path of the folder where the model output, including the summary and predictions, are saved"
	outputpath::TS=""
	"event on each trial aligned to which spikes are counted"
	reference_event::TS="cpoke_in"
	"time, in second, after which spikes are included on each trial, aligned to the reference event on that trial."
	time_in_trial_begin_s::TF=-2.0
	"time, in second, before which spikes are included on each trial, aligned to the reference event on that trial."
	time_in_trial_end_s::TF=2.0; @assert time_in_trial_end_s > time_in_trial_begin_s
	"absolute path to a MAT file containing the indices of the trials to be used"
	trial_indices_path::TS=""
	"event on each trial after which spikes are not counted"
	trim_after_event::TS=""
end

"""
    Trial

Sensory stimuli, behavior, and spike train on each trial
"""
@with_kw struct Trial{	TB<:Bool,
						TF<:AbstractFloat,
						TI<:Integer,
						TVI<:Vector{<:Integer},
						TARF<:AbstractRange{<:AbstractFloat},
						VVF<:Vector{<:Vector{<:AbstractFloat}}}
	"behavioral choice"
    choice::TB
	"source of each click: `0` indicates a left click, `1` a right click, and `2` a stereoclick"
	clicks_source::TVI
	"time step of each click"
	clicks_timestep::TVI
	"first complete trial in the session"
	first_reference_time_s::TF
	"log of the ratio of the generative right and left click rate"
	γ::TF
	"last complete trial in the session"
	last_reference_time_s::TF
	"time step of leaving the center port"
	movement_timestep::TI
	"pose measurements"
	pose::VVF
	"time of the reference event"
	reference_time_s::TF
	"time step of entering the side port"
	response_timestep::TI
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
	"number of time steps before the trial used for estimating the post-spike filter"
	Tpre::TI = length(ypre)
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
	baseline::UI
	click::UI
	stereoclick::UI
	leftclick::UI
	leftmovement::UI
	leftresponse::UI
	movement::UI
	pose::UI
	postspike::UI
	response::UI
	rightclick::UI
	rightmovement::UI
	rightresponse::UI
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
	"inferred baseline firing rate on each trial"
	baseline::VF
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
	"concatenated weights"
	𝐰_baseline::VF
	"precision parameter of the Gaussian prior on weights"
	a::VF=[options.precision]
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
@with_kw struct EvidenceOptimization{VR<:Vector{<:Real}, VVR<:Vector{<:Vector{<:Real}}, VMF<:Vector{<:Matrix{<:AbstractFloat}}}
	"precision used for each iteration of MAP optimization"
	a::VR
	"approximate log-evidence evaluated at the end of each MAP optimization"
	𝐸::VR
	"a matrix of the second-order partial derivatives of the log-likelihood with respect to the model parameters"
	hessian_loglikelihood::VMF
	"Euclidean norm of gradient of the log posterior at the end of each MAP optimization"
	MAP_g_residual::VR
	"MAP parameters"
	𝐰::VVR
end

"""
	PerieventTimeHistogram

The mean across trials of a single condition (e.g. trials that ended with a left choice) of the filtered spike train of one neuron, and the estimated 95% confidence interval of the trial mean
"""
@with_kw struct PerieventTimeHistogram{S<:String, VF<:Vector{<:AbstractFloat}}
	"condition"
	condition::S
	"estimate of the peri-event time histogram based on observed spike trains"
	observed::VF
	"estimate of the peri-event time histogram based on simulated spike trains"
	predicted::VF
	"event on each trial to which the histogram is aligned"
	reference_event::S
	"time, in seconds"
	timesteps_s::VF
end

"""
	PETHSet

A set of peri-event time histogram of one neuron.
"""
@with_kw struct PETHSet{PETH<:PerieventTimeHistogram}
	"average across trials that ended in a left choice, including both correct and incorrect trials"
	leftchoice::PETH
	"average across trials on which the aggregate evidence favored left and the reward is baited on the left"
	leftevidence::PETH
	"average across trials on which the animal made a left choice and the evidence was strongly leftward. Therefore, only correct trials are included. Evidence strength is defined by the generative log-ratio of click rates: `γ` ≡ log(right click rate) - log(left click rate). Strong left evidence evidence include trials for which γ < -2.25"
	leftchoice_strong_leftevidence::PETH
	"average across trials on which the animal made a left choice and the generatative `γ` < -2.25 && `γ` < 0"
	leftchoice_weak_leftevidence::PETH
	rightchoice::PETH
	rightevidence::PETH
	rightchoice_strong_rightevidence::PETH
	rightchoice_weak_rightevidence::PETH
	"average across all trials"
	unconditioned::PETH
end

"""
	Kernel

Linear filter of signal from task events or past spiking
"""
@with_kw struct Kernel{VF<:Vector{<:AbstractFloat}, S<:String}
	"to identify the set of basis functions"
	basisname::S
	"the filter"
	h::VF
	"to identify the input to which this filter is associated"
	inputname::S
	"the time of each value of the filter `h` aligned to the time of the input"
	timesteps_s::VF
end

"""
	Characterization
"""
@with_kw struct Characterization{VI<:Vector{<:Integer}, VVI<:Vector{<:Vector{<:Integer}}, VVR<:Vector{<:Vector{<:Real}}}
	"log-likelihood (nats) on each time step on each trial"
	LL::VVR
	"total pre-activation input that is external to the neuron (i.e., no spike history)"
	externalinput::VVR
	"moment-to-moment firing rate inferred as the expectation across simulations"
	inferredrate::VVR
	"autocorrelation function"
	autocorrelation::VVR
	"moment-to-moment spike counts"
	observed_spiketrains::VVI
	"information for identifying the data"
	trialindices::VI
end

"""
	CVIndices

Indices of trials and timesteps used for training and testing
"""
@with_kw struct CVIndices{VI<:Vector{<:Integer}}
	"indices of the trials used for training or testing"
	testingtrials::VI
	trainingtrials::VI
end

"""
	CVResults

Results of cross-validation
"""
@with_kw struct CVResults{C<:Characterization, VC<:Vector{<:CVIndices}, VS<:Vector{<:Model}, VPETH<:Vector{<:PerieventTimeHistogram}}
	"a composite containing quantities that are computed out-of-sample and used to characterize the model`"
	characterization::C
	"cvindices[k] indexes the trials and timesteps used for training and testing in the k-th resampling"
	cvindices::VC
	"peri-event time histograms"
	peths::VPETH
	"training models"
	trainingmodels::VS
end
