"""
RETURN an instance of a type converted into an instance of a `Dict`
"""
dictionary(x) = Dict((String(fieldname)=>getfield(x,fieldname) for fieldname in fieldnames(typeof(x)))...)

"""
"""
function save(characterization::Characterization, eo::EvidenceOptimization, model::Model)
end
