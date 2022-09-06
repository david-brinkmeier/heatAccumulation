classdef optsheatsource
    % This object contains special settings for setup of heat sources
    % At the moment this only applies to generation of line heat sources from
    % superposition of point sources.
    % These settings are parsed for example to get2D_Distributed_Heatkernels.m
    % for calculation of heat kernels
    % https://de.mathworks.com/help/matlab/matlab_oop/validate-property-values.html
    % https://de.mathworks.com/help/matlab/matlab_oop/property-validator-functions.html
    
    properties
        % all lengths in SI-units
        name                string          % identifier
        sz_heat(:,1)        cell            % beam radius in SI units (0 = point source, != 0 Gaussian source) or array with corresponding radius to lheat
        x_pos(1,1)          double          % lateral X position where evaluation occurs
        y_pos(1,1)          double          % lateral Y Position where evaluation occurs
        z_pos(1,1)          double          % 0 = surface, positive values = inside the borehole
                                            % only scalar allowed, evaluate different z-positions using loops
        l_heat(:,1)         cell            % scalar or vector or cell array of scalar or cell array of vector
                                            % scalar: heat sources will be deployed from 0 to scalar, if vector then vector describes position of sources
        q_vect(:,1)         cell            % relative energy distribution along z: for details see get2D_Distributed_Heatkernels.m
                                            % may be string/char, 1D vector or 2D Matrix of length (length(lheat),arbitrary)
        n_heatsources(1,1)  double          % used if no variable spacing using delta_l_heat is used
        delta_l_heat(1,1)   double          % specifies maximum spacing between point sources
    end
    properties (SetAccess = protected)
        % placeholder
    end
    properties (Dependent) %(Dependent, SetAccess = private)
        % placeholder
    end
    
    methods
        % constructor
        function obj = optsheatsource(varargin)
            % if a single string is parsed we compare and initialize
            if ~isempty(varargin) && numel(varargin) == 1
                obj.name = varargin{1}; % set Name
            else
                obj.name = 'point';     % default to points
            end
            obj.sz_heat = {0};          % size of the heat source (for gaussian source equal to beam radius)
            obj.x_pos = 0;              % default position: on axis
            obj.y_pos = 0;              % default position: on axis
            obj.z_pos = 0;              % default position: surface
            obj.l_heat = {0};           % this must be set by user
            obj.q_vect = {'constant'};  % default to constant energy distribution
            obj.n_heatsources = 100;    % pulse length in s
            obj.delta_l_heat = nan;      % maximum distance between heatsources
        end
        
        function obj = set.name(obj,value)
            allowed = {'point','gauss','ring'};
            value = lower(char(value));
            if ismember(value,allowed)
                obj.name = value;
            else
                error('Heatsource "%s" is not defined for optsheatsource',value);
            end
        end
        
        function obj = set.n_heatsources(obj,value)
            obj.n_heatsources = value;
            % if user sets n_heatsources assume delta_l_heat should not be used and clear it
            obj = triggerReset(obj,'delta_l_heat');
        end
        
        function value = get.sz_heat(obj)
            if strcmpi(obj.name,'gauss')
                % bauer / michalowski gaussian heat source produces
                % singularity, cannot allow zeros, only eps.
                value = obj.fixzeros(obj.sz_heat);
            else
                value = obj.sz_heat;
            end
        end
    end

%% static / private methods

    methods (Access = private)
        function obj = triggerReset(obj,valname)
            % resets value with name valname
            obj.(valname) = nan;
        end
    end
    
    methods (Static)
        function x = fixzeros(x)
            % recursively replace zeros by eps for numerics / cells
            if isnumeric(x)
                x(x==0) = eps;
            elseif iscell(x)
                x = cellfun(@optsheatsource.fixzeros,x,'un',0);
            end
        end        
    end
end
