function [ref, varargout] = denormalize_data(ref, factor, varargin)
    ref = ref * factor; 
    
    for i = 1:length(varargin)
        varargout{i} = varargin{i} * factor; 
    end

end