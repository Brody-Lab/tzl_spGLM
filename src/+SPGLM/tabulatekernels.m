function T = tabulatekernels(kernels)
%{
RETURN a table of the kernels based on the results from the SPGLM Julia module

ARGUMENT
-kernels: a cell vector whose each element is a struct containing information for a convolution
filter for an input
%}
validateattributes(kernels, {'cell'}, {'vector'})
field_names = sort(string(fieldnames(kernels{1}))');
T = struct;
for name = field_names
    for i = 1:size(kernels,1)
        T.(name){i,1} = kernels{i}.(name);
    end
    if iscellstr(T.(name))
        T.(name) = string(T.(name));
    end
end
T.inputname(T.inputname == "time_in_trial") = "fixation";

T = struct2table(T);