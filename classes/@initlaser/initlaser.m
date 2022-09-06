classdef initlaser
    % Generate structure for (current) laser source properties
    
    properties
        name        string % name / identifier of the laser source
        Ep          double % pulse energy in Joule
        frep        double % pulse repetition rate in Hz
        w0          double % gaussian beam radius in m
        tp          double % pulse length in s
    end
    properties (Dependent) %(Dependent, SetAccess = private)
        tpp         double % pulse to pulse length
    end
    
    methods
        % constructor
        function obj = initlaser(varargin)
            % if a single string is parsed we compare and initialize a
            % laser (if it is known)
            if ~isempty(varargin) && numel(varargin) == 1
                if strcmpi('trumicro',varargin{:})
                    obj.name = 'trumicro 1030nm'; % set Name
                    obj.Ep = 0;             % pulse energy in Joule
                    obj.frep = 3e5;         % pulse repetition rate in Hz
                    obj.w0 = 25e-6;         % gaussian beam radius in m
                    obj.tp = 7e-12;         % pulse length in s
                end
            else
                obj.name = 'unknown';
            end
        end
        
        % calculate values of dependent properties
        function value = get.tpp(obj)
            % calculate pulse-to-pulse length from reprate
            value = 1/obj.frep;
        end
        % error message when user attempts to set dependent property tpp
        function obj = set.tpp(obj,~)
            warning('tpp property is automatically calculated from pulse repetition rate\n');
            fprintf('%s%d\n','Pulse to pulse length tpp in µs is: ',obj.tpp*1e6)
        end
    end
end