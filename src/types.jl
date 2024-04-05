"""
    Options

Model settings
"""
@with_kw struct Options{TF<:AbstractFloat, TS<:String}
	"full path of the data"
    datapath::TS=""
	"duration of each timestep in seconds"
    Δt::TF=0.01
	"absolute path of the folder where the model output, including the summary and predictions, are saved"
	outputpath::TS=""
end
