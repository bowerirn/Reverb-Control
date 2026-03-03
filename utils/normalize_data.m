function [ref, factor, varargout] = normalize_data(ref, varargin)
    factor = max(abs(ref));
    ref = ref / factor; % Normalize the wav data
    
    for i = 1:length(varargin)
        varargout{i} = varargin{i} / factor; % Optional output based on varargin
    end

end