"""
RETURN an instance of a type converted into an instance of a `Dict`
"""
dictionary(x) = Dict((String(fieldname)=>getfield(x,fieldname) for fieldname in fieldnames(typeof(x)))...)
function dictionary(eo::EvidenceOptimization)
    Dict("a"=>eo.a,
        "E"=>eo.𝐸,
        "MAP_g_residual"=>eo.MAP_g_residual,
        "w"=>eo.𝐰)
end

function dictionary(characterization::Characterization)
    Dict("LL"=>characterization.LL,
          "inferredrate"=>characterization.inferredrate,
          "peths"=>map(dictionary, characterization.peths))
end

function dictionary(basisset::BasisFunctionSet)
    Dict("name"=>String(basisset.name),
          "timesteps_s"=>collect(basisset.timesteps_s),
          "Phi"=>basisset.Φ)
end
function dictionary(model::Model)
    Dict("a"=>model.a,
          "options"=>dictionary(model.options),
          "basissets"=>map(dictionary, model.basissets),
          "weightindices"=>Dict((String(fieldname)=>collect(getfield(model.weightindices,fieldname)) for fieldname in fieldnames(typeof(model.weightindices)))...),
          "w"=>model.𝐰)
end
"""
WRITE a MATLAB ".mat" file
"""
save(model::Model) = matwrite(joinpath(model.options.outputpath, "model.mat"), dictionary(model))
save(characterization::Characterization, outputpath::String) = matwrite(joinpath(outputpath, "characterization.mat"), dictionary(characterization))
save(eo::EvidenceOptimization, outputpath::String) = matwrite(joinpath(outputpath, "evidenceoptimization.mat"), dictionary(eo))
