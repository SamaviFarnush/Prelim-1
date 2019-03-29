
function [parameter_array] = getparam(end_index,parameter,delta)
    
    parameter_array = zeros(end_index,end_index);
    for index = 1:end_index
        load_parameter = parameter;
        new_parameter = load_parameter(index)*(1+delta);
        load_parameter(index) = new_parameter;
        parameter_array(index,:) = load_parameter;
    end
end