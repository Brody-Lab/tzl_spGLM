function T = tabulatepeths(peths)
%{
tabulate peri-event time histograms generated by the Julia module SPGLM

ARGUMENT
-`peths`: a cell vector

OUTPUT
-a table
%}
validateattributes(peths, {'cell'}, {'vector'})
field_names = sort(string(fieldnames(peths{1}))');
T = struct;
for name = field_names
    for i = 1:size(peths,1)
        T.(name){i,1} = peths{i}.(name);
    end
    if iscellstr(T.(name))
        T.(name) = string(T.(name));
    end
end

T = struct2table(T);