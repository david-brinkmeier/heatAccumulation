classdef simopts
    % Generate structure to store relevant material constants for a material
    % types of properties: https://www.mathworks.com/help/matlab/matlab_oop/mutable-and-immutable-properties.html
    % https://www.mathworks.com/help/matlab/matlab_oop/how-to-use-properties.html#brh76wy-1_2
    % todo: consider generating q_distribFull implicity / serve through public method
    
    properties
        UUID                    char            % identifier
        oversampling(1,2)       cell            % enable if you want to resolve Tempearature in between pulses, recommed integer >= 5
        TotalNoPulses(1,1)      double          % for numerical efficiency this is ideally (2^n)/2 where n = positive integer, e.g. 1,2,3...
        bins(1,2)               cell            % [2^n,'static'] or 'dynamic': Binning to reduce numeric complexity. Chunksize = TotalNoPulses/bins (if 'static')
        PulseRange(1,2)         double          % vector [firstpulse, lastpulse] start pulsing at firstpulse, stop at lastpulse [but temporally evaluated until TotalNoPulses]
        EnergyRamp(1,:)         double          % vector of length(1,TotalNoPulses) which contains incident energy scalars (for energy ramps)
        Depths_Eval(1,2)        cell            % if applicable (i.e. currently only for "3DLine" enter vector of to-be-evaluated z-positions
                                                % or cell such as {'string',number}, ex. {'determine',50} evaluates 50 linearly spaced positions from 0 to max drill depth
        radial_position(1,:)    double          % (1,2) Vector of X and Y Position where Evaluation occurs in SI-units
        CorrfactAbs(1,:)        double          % Energy correction factor of length (TotalNoPulses,1) applied as energy_in = corrfact(i)*mat.eta_abs*laser.Ep
        Heatsource2Dlen(:,1)    cell            % Cell (TotalNoPulses,1) containing length of the heat source at/after each Pulse or vector
        Heatsource2Dsize(:,1)   cell            % Cell (TotalNoPulses,1) containing radial size of heat source at/after each Pulse or vector (applies to 3DLines with gaussian source)
        Q_distrib_2D(:,1)       cell            % cell containing vector (1,:), string / char (1,1) or 2D Matrix (TotalNoPulses,:) or array of cells (TotalNoPulses,:)) 
                                                % for despription of energy distribution along 2D heat source from 3D point sources
        Heatsource              char            % Select type of heatsource ("3Dgauss" = Gaussian heatsource, "3D" = 3d point source, "2D" = variable length line sources, "3DLine" line source from point sources)
        HeatsourceSettings      optsheatsource  % containts heat source specific settings (currently only for 3D distributed point sources)
        precision(1,2)          cell            % {computation,storage}; fp64-fp32-fp16; for calculation fp32 is min, storage is applied at last step to avoid repeated roundoff errors
        isdrilling(1,2)         double          % [flag,ablcrater]; drillmodel must set this flag and provide ablation radius
                                                % isdrilling(1) is flag, isdrilling(2) hold ablation crater radius
                                                % proper implementation when better drillmodels are available that can handle varying pulse energies etc
        debug                   logical         % adds some plots for drill depth progress and heat accumulation
    end
    properties (SetAccess = protected)
        PulsesPerBin                double      % number of pulses per bin; Scalar OR Matrix containing indices / pulses per bin derived from Heatsource2Dlen
        z_eval_depths               double = 0  % calculated dependent on input
        timesteps(1,3)              cell        % containts name (constant/variable) / timesteps vector normal / oversampled (times vector from t0 to tEnd)
        timestepping(1,2)           cell        % actual timesteps between consecutive pulses normal / oversampled
        reprates(1,2)               cell        % applied reprates / applied reprates corresponding to oversampled results
        corrfactBin(1,:)            double      % binned energy correction factor (1,sim.bins) calculated from CorrfactAbs for convolution calc
        corrfactFull(1,2)           cell        % extrapolated "true" corrfact from corrfactBin for explicit calc
        heatsource2DlenBin(:,1)     cell        % binned heat source 2D length (1,sim.bins) calculated from Heatsource2Dlen for convolution calc
        heatsource2DlenFull(1,2)    cell        % extrapolated "true" heat source 2D length from heatsource2DlenBin for explicit calc
        heatsource2DsizeBin(:,1)    cell        % binned sizes of the heat source (currently only applies to gaussian heat sources)
        heatsource2DsizeFull(1,2)   cell        % extrapolated "true" sizes of source from heatsource2DsizeBin for explicit calc
        q_distribBin(:,1)           cell        % calculated from Q_distrib_2D: binned energy distribution along 2D heat source from 3D point sources
        q_distribFull(1,2)          cell        % extrapolated "true" energy distribution for explicit calc
        reprate_max(1,1)            double      % grab maximum reprate from chirp vector or (if constant) as provided by laser.frep
    end
    properties (Dependent)
        EnergyRampBin(1,:)    double            % LEGACY EnergyRampBin vector of length [generates pulsematrix (kron) matrix with heatinputs_per_bin normal / oversampled]
        HeatInputsPerBin(1,2) cell              % Cell array containting normal and oversampled HeatInputsPerBin
        EnergyRampFull(1,2)   cell              % [Scaled Heatinput Vector] derived from EnergyRampBin: We need to use this - Drill depth and binned heat accumulation MUST use the same input
                                                % Oversampled version is derived from EnergyRampFull: Required if we calculate oversampled heat accumulation with energy ramp
        HeatInputs(1,2)       cell              % Cell array containing normal oversampled HeatInputs vector
    end
    
    methods
        % Constructor
        function obj = simopts(varargin)
            if ~isempty(varargin) && numel(varargin) == 1
                if strcmp('init',varargin{:})
                    obj.UUID               = obj.generateID(5);
                    obj.oversampling       = {1,1};
                    obj.TotalNoPulses      = 2^10;
                    obj.PulseRange         = [1, 2^10];
                    obj.EnergyRamp         = ones(1,2^10);
                    obj.Q_distrib_2D       = {'constant'};
                    obj.Heatsource2Dsize   = {0};
                    obj.radial_position    = 0;
                    obj.Depths_Eval        = {0e-6,[]};
                    obj.Heatsource         = 'point';
                    obj.HeatsourceSettings = optsheatsource(obj.Heatsource);
                    obj.precision          = {'fp32','fp32'};
                    obj.isdrilling         = [0,0];                    
                    obj.debug              = true;
                    obj.useallcores()
                end
            else
                warning('initialize simopts using "init" flag')
            end
        end  
        
        %% GET METHODS
        
        % Calculates dependent property bins
        function value = get.bins(obj)
            % bins specifies the number of "chunks" that the heatkernels /
            % energy distribution / size of heat sources (if not point source)
            % will be divided into. if bins << totalnopulses then
            % 'variable' bin length helps to suppress aberrations longer by
            % smarter distribution of pulses per bin
            %
            % if n_pulses < bins{1} then we set bins = n_pulses
            % if n_pulses = bins{1} this means one heat kernel per pulse
            % if n_pulses > bins{1} this means one bin contains PulsesPerBin = n_pulses/bins pulses
            value = obj.bins;
            % we can't have more bins than pulses
            if value{1} > obj.TotalNoPulses
                value{1} = obj.TotalNoPulses;
                warning('Bins cannot exceed total number of pulses, fixed.')
            end
            % we can't have anything else than a power of two
            % this is only triggered if the case above has not been triggered
            if ~(2^nextpow2(value{1}) == value{1})
                value{1} = 2^nextpow2(obj.bins{1});
                warning(strcat('bins{1} is not 2^n, bins has been set to',32,num2str(value{1})))
            end
        end
                
        function value = get.PulseRange(obj)
            % sanity check on PulseRange limits
            range = obj.PulseRange;
            if range(1) < 1
                warning('Pulserange cannot have negative values, Pulserange(1) set to 1')
                range(1) = 1;
            end
            if range(2) > obj.TotalNoPulses
                warning('Pulserange cannot exceed TotalNoPulses, Pulserange(2) set to TotalNoPulses')
                range(2) = obj.TotalNoPulses;
            end
            value = range;
        end
        
        function value = get.EnergyRamp(obj)
            % if nothing 
            value = obj.EnergyRamp;
            if length(value) ~= obj.TotalNoPulses
                warning('Energy ramp is set to default: Disabled / constant power. Input must be vector (1,TotalNoPulses)')
                value = ones(1,obj.TotalNoPulses);
            end
            % all values after PulseRange2 are laser OFF
            if obj.PulseRange(2) > 0
                value(1,obj.PulseRange(2):end) = 0;
            else
                warning('PulseRange cannot have negative values -> ignored')
            end
        end
        
        function value = get.EnergyRampBin(obj)
            % the binned vector is simply an undersampled version of the
            % user-provided original (if binning is used)
            % only used with legacy / kron pulsematrix / static binning
            if strcmp(obj.bins{2},'static')
                value = obj.EnergyRamp(1:obj.PulsesPerBin:end);
            else
                value = NaN; % only defined / applicable for 'static'
            end
        end
        
        function value = get.EnergyRampFull(obj)
            value{1} = obj.EnergyRamp;
            % Oversampled variant / Required to calculate ramped heat input vector
            value{2} = repelem(value{1},obj.oversampling{2});
        end
        
        function value = get.HeatInputsPerBin(obj)
            % Required to generate (binned) sparse Pulsematrix
            % First pulse is always on with binning
            % only used with legacy / kron pulsematrix / static binning
            if strcmp(obj.bins{2},'static')
                value{1} = ones(1,obj.PulsesPerBin);
                value{2} = zeros(1,obj.oversampling{2}*obj.PulsesPerBin);
                value{2}(obj.oversampling{2}:obj.oversampling{2}:obj.PulsesPerBin*obj.oversampling{2}) = 1;
            else
                value = {NaN,NaN}; % only defined / applicable for 'static'
            end
        end
        
        function value = get.HeatInputs(obj)
            % used whenever binning is not applied
            % these vectors represent boolean heat inputs [normal / oversampled]
            value{1} = ones(1,obj.TotalNoPulses);
            % oversampled has "cooling periods" / timesteps evaluated
            % between heat inputs
            value{2} = zeros(1,obj.oversampling{2}*obj.TotalNoPulses);
            value{2}(obj.oversampling{2}:obj.oversampling{2}:obj.TotalNoPulses*obj.oversampling{2}) = 1;
        end
        
        function value = get.HeatsourceSettings(obj)
            % generate nested class HeatsourceSettings based on whether
            % an appropriately generated object currently exists
            if ~strcmpi(obj.Heatsource,obj.HeatsourceSettings.name)
                oldvals = {obj.HeatsourceSettings.n_heatsources,obj.HeatsourceSettings.delta_l_heat};
                value = optsheatsource(obj.Heatsource); % init with name
                % write old vals
                value.n_heatsources = oldvals{1};
                value.delta_l_heat = oldvals{2}; 
            else
                value = obj.HeatsourceSettings;
            end
        end
        
        function value = get.CorrfactAbs(obj)
            % CorrfactAbs is meant to account for increased absorption, i.e.
            % through Gouffé or Raytracing. For each pulse (geometry) one
            % scalar should be provided, i.e. CorrfactAbs is a vector
            % (1,TotalNoPulses) or (TotalNoPulses,1)
            if (length(obj.CorrfactAbs) ~= obj.TotalNoPulses)
                % not set? make it all ones!
                value = ones(1,obj.TotalNoPulses);
            else
                value = obj.CorrfactAbs;
            end
        end
        
        function value = get.corrfactBin(obj)
            % take mean of every PulsesPerBin values
            % if it is not empty and it is of correct length then do generate bin vector
            % only used with legacy / kron pulsematrix / static binning
            if strcmp(obj.bins{2},'static')
                value = mean(reshape(obj.CorrfactAbs,[obj.PulsesPerBin, obj.bins{1}]),1);
            else
                value = NaN; % only defined / applicable for 'static'
            end
        end
        
        function value = get.corrfactFull(obj)
            % generate corrfactFull for internal usage
            value = {repelem(obj.CorrfactAbs,obj.oversampling{1}),...
                     repelem(obj.CorrfactAbs,obj.oversampling{2})};
        end
        
        function value = get.heatsource2DlenBin(obj)
            % if we call heatsource2DlenBin then we need a Heatsource2Dlen
            % in any case we need enough cells so that we heave one length
            % corresponding to each pulse
            % alternatively allowed are vectors for each pulse where each
            % vector represents depths where heat sources will be placed
            if ~isequal(size(obj.Heatsource2Dlen),[obj.TotalNoPulses 1])
                error('Heatsource2Dlen must be cell (TotalNoPulses,1)')
            end
            % take mean of every PulsesPerBin values
            if max(cellfun(@(x) size(x,2), obj.Heatsource2Dlen)) == 1
                % then each cell contains only a scalar. convert to matrix,
                % take mean of every obj.PulsesPerBin, convert to cell
                if strcmp(obj.bins{2},'variable')
                    % then obj.PulsesPerBin is a matrix and Pulses per Bin
                    % are variable according to indices stored in PulsesPerBin
                    value = obj.Heatsource2Dlen(obj.PulsesPerBin(:,3)); % just extract corresponding to mean idx
                    % plot(obj.PulsesPerBin(:,3),cell2mat(value)), hold on, plot(cell2mat(obj.Heatsource2Dlen)) % check
                    %
                    % below alternative more elaborate variant (consider implementing later)
                    % hslenghts = cell2mat(obj.Heatsource2Dlen); % store heatsource lengths
                    % value = cell(obj.bins{1},1); % init
                    % for i = 1:obj.bins{1}; value(i) = {mean(hslenghts(obj.PulsesPerBin(i,1):obj.PulsesPerBin(i,2)))}; end
                else
                    value = num2cell(mean(reshape(cell2mat(obj.Heatsource2Dlen),[obj.PulsesPerBin, obj.bins{1}]).',2),2);
                end
            else
                % then Heatsource2Dlen must be vectors and each row describes
                % the positions of the heat sources for that pulse
                % value = obj.Heatsource2Dlen(1:obj.PulsesPerBin:end);
                value = [obj.Heatsource2Dlen(floor(1+obj.PulsesPerBin/2):obj.PulsesPerBin:floor(end-obj.PulsesPerBin/2)); obj.Heatsource2Dlen(end)];
                if ~isequal(size(value),[obj.bins{1} 1]), error('generated heatsource2DlenBin has wrong size'); end
            end
        end
        
        function value = get.heatsource2DlenFull(obj)
            % if Heatsource2Dlen is bad defined then heatsource2DlenBin will throw error
            % since heatsource2DlenFull is derived from bin vector
            if strcmp(obj.bins{2},'variable')
                value = {cell(obj.TotalNoPulses*obj.oversampling{1},1),...
                         cell(obj.TotalNoPulses*obj.oversampling{2},1)}; % init
                hslenbin = obj.heatsource2DlenBin; % avoid repeated calls to dependent var
                for i = 1:obj.bins{1}
                    if obj.PulsesPerBin(i,4)
                        % then repeat the elements according to the indices stored in PulsesPerBin(i,4) which means
                        % indices obj.PulsesPerBin(i,1) through obj.PulsesPerBin(i,2) (obj.PulsesPerBin(i,4) pulses times oversampling 1)
                        value{1}(obj.PulsesPerBin(i,1):obj.PulsesPerBin(i,2),1) = repelem(hslenbin(i),obj.PulsesPerBin(i,4)*obj.oversampling{1}).';
                        % indices obj.PulsesPerBin(i,5) through obj.PulsesPerBin(i,6) (obj.PulsesPerBin(i,4) pulses times oversampling 2)
                        value{2}(obj.PulsesPerBin(i,5):obj.PulsesPerBin(i,6),1) = repelem(hslenbin(i),obj.PulsesPerBin(i,4)*obj.oversampling{2}).';
                    end
                end
            else
                value = {repelem(obj.heatsource2DlenBin,obj.PulsesPerBin*obj.oversampling{1}),...
                    repelem(obj.heatsource2DlenBin,obj.PulsesPerBin*obj.oversampling{2})};
            end
        end
        
        function value = get.Heatsource2Dsize(obj)
            % if Heatsource2Dsize is a scalar we repeat if we
            % obj.TotalNoPulses is provided
            if isscalar(obj.Heatsource2Dsize)
                value = repmat(obj.Heatsource2Dsize,[obj.TotalNoPulses 1]);
            elseif isequal(size(obj.Heatsource2Dsize),[obj.TotalNoPulses 1])
                % else sizes are provided for each pulse
                % these sizes may be scalars (constant size for all axial
                % positions for that pulse OR vectors where the heat source
                % size vector describes the size of the heatsource from 0
                % (top of workpiece) to Heatsource2DLen(end,nPulse)
                value = obj.Heatsource2Dsize;
            else
                error('Heatsource2Dsize must either be cell (1,1) or cell(TotalNoPulses,1)')
            end
        end

        function value = get.heatsource2DsizeBin(obj)
            % if we call heatsource2DlenBin then we need a Heatsource2Dlen
            % in any case we need enough cells so that we heave one length
            % corresponding to each pulse
            % alternatively allowed are vectors for each pulse where each
            % vector represents depths where heat sources will be placed
            if ~isequal(size(obj.Heatsource2Dsize),[obj.TotalNoPulses 1])
                error('Heatsource2Dlen must be cell (TotalNoPulses,1)')
            end
            % take mean of every PulsesPerBin values
            if max(cellfun(@(x) size(x,2), obj.Heatsource2Dsize)) == 1
                % then each cell contains only a scalar. convert to matrix,
                % take mean of every obj.PulsesPerBin, convert to cell
                if strcmp(obj.bins{2},'variable')
                    value = obj.Heatsource2Dsize(obj.PulsesPerBin(:,3)); % just extract corresponding to mean idx
                else
                    value = num2cell(mean(reshape(cell2mat(obj.Heatsource2Dsize),[obj.PulsesPerBin, obj.bins{1}]).',2),2);
                end
            else
                % then Heatsource2Dsize must be vectors and each row describes
                % the sizes of the heat sources corresponding to the heatsource2Dlen
                if strcmp(obj.bins{2},'variable')
                    value = obj.Heatsource2Dsize(obj.PulsesPerBin(:,3)); % just extract corresponding to mean idx
                else
                    value = [obj.Heatsource2Dsize(floor(1+obj.PulsesPerBin/2):obj.PulsesPerBin:floor(end-obj.PulsesPerBin/2)); obj.Heatsource2Dsize(end)];
                end
            end
            if ~isequal(size(value),[obj.bins{1} 1]), error('generated heatsource2DsizeBin has wrong size'); end
        end
        
        function value = get.heatsource2DsizeFull(obj)
            % if Heatsource2Dlen is bad defined then heatsource2DlenBin will throw error
            % since heatsource2DlenFull is derived from bin vector
            if strcmp(obj.bins{2},'variable')
                value = {cell(obj.TotalNoPulses*obj.oversampling{1},1),...
                         cell(obj.TotalNoPulses*obj.oversampling{2},1)}; % init
                hsszbin = obj.heatsource2DsizeBin; % avoid repeated calls to dependent var
                for i = 1:obj.bins{1}
                    if obj.PulsesPerBin(i,4)
                        % then repeat the elements according to the indices stored in PulsesPerBin(i,4) which means
                        % indices obj.PulsesPerBin(i,1) through obj.PulsesPerBin(i,2) (obj.PulsesPerBin(i,4) pulses times oversampling 1)
                        value{1}(obj.PulsesPerBin(i,1):obj.PulsesPerBin(i,2),1) = repelem(hsszbin(i),obj.PulsesPerBin(i,4)*obj.oversampling{1}).';
                        % indices obj.PulsesPerBin(i,5) through obj.PulsesPerBin(i,6) (obj.PulsesPerBin(i,4) pulses times oversampling 2)
                        value{2}(obj.PulsesPerBin(i,5):obj.PulsesPerBin(i,6),1) = repelem(hsszbin(i),obj.PulsesPerBin(i,4)*obj.oversampling{2}).';
                    end
                end
            else
                value = {repelem(obj.heatsource2DsizeBin,obj.PulsesPerBin*obj.oversampling{1}),...
                    repelem(obj.heatsource2DsizeBin,obj.PulsesPerBin*obj.oversampling{2})};
            end
        end
        
        function value = get.q_distribBin(obj)
            if isequal(size(obj.Q_distrib_2D),[1 1])
                % then it is [1,1] cell
                % it contains a valid matrix if
                if size(obj.Q_distrib_2D{:},1) == obj.TotalNoPulses
                    % in which case the matrix must be binned so that
                    value = obj.binqvectmatrix(cell2mat(obj.Q_distrib_2D),obj.bins,obj.TotalNoPulses,obj.PulsesPerBin);
                elseif size(obj.Q_distrib_2D{:},1) == 1
                    % if it is no matrix but 1x1 cell then should be 1D vector or string
                    % interpretation of strings is handled by heatkernel calculation / optsheatsource class
                    value = repmat(obj.Q_distrib_2D,obj.bins{1},1);
                else
                    error('Q distribution must either be char / 1D or cell(TotalNoPulses,1) or matrix (TotalNoPulses,~)')
                end
            elseif isequal(size(obj.Q_distrib_2D),[obj.TotalNoPulses 1])
                % then it is array of cells (TotalNoPulses,1)
                if strcmp(obj.bins{2},'variable')
                    value = obj.Q_distrib_2D(obj.PulsesPerBin(:,3),:); % just extract corresponding to mean idx
                else
                    value = [obj.Q_distrib_2D(floor(1+obj.PulsesPerBin/2):obj.PulsesPerBin:floor(end-obj.PulsesPerBin/2)); obj.Q_distrib_2D(end)];
                end
            else
                error('Q_distrib_2D must either be cell (1,1) or cell(TotalNoPulses,1)')
            end
        end
        
        function value = get.q_distribFull(obj)
            if isequal(size(obj.Q_distrib_2D),[1 1])
                % then it is [1,1] cell
                % it contains a valid matrix if
                if size(obj.Q_distrib_2D{:},1) == obj.TotalNoPulses
                    % in which case the matrix must be binned so that
                    value = obj.unbinqvectmatrix(cell2mat(obj.q_distribBin),obj.bins,obj.TotalNoPulses,obj.PulsesPerBin,obj.oversampling); % for explicit
                elseif size(obj.Q_distrib_2D{:},1) == 1
                    % if it is no matrix then it is 1D vector or string
                    value = {repmat(obj.Q_distrib_2D,obj.TotalNoPulses*obj.oversampling{1},1),...
                             repmat(obj.Q_distrib_2D,obj.TotalNoPulses*obj.oversampling{2},1)};
                else
                    error('Q distribution must either be char / 1D or cell(TotalNoPulses,1) or matrix (TotalNoPulses,~)')
                end
            elseif isequal(size(obj.Q_distrib_2D),[obj.TotalNoPulses 1])
                % else it should be array of cells (TotalNoPulses,1)
                if strcmp(obj.bins{2},'variable')
                    value = {cell(obj.TotalNoPulses*obj.oversampling{1},1),...
                             cell(obj.TotalNoPulses*obj.oversampling{2},1)}; % init
                    qdistribin = obj.q_distribBin; % avoid repeated calls to dependent var
                    for i = 1:obj.bins{1}
                        if obj.PulsesPerBin(i,4)
                            % then repeat the elements according to the indices stored in PulsesPerBin(i,4) which means
                            % indices obj.PulsesPerBin(i,1) through obj.PulsesPerBin(i,2) (obj.PulsesPerBin(i,4) pulses times oversampling 1)
                            value{1}(obj.PulsesPerBin(i,1):obj.PulsesPerBin(i,2),1) = repelem(qdistribin(i),obj.PulsesPerBin(i,4)*obj.oversampling{1}).';
                            % indices obj.PulsesPerBin(i,5) through obj.PulsesPerBin(i,6) (obj.PulsesPerBin(i,4) pulses times oversampling 2)
                            value{2}(obj.PulsesPerBin(i,5):obj.PulsesPerBin(i,6),1) = repelem(qdistribin(i),obj.PulsesPerBin(i,4)*obj.oversampling{2}).';
                        end
                    end
                else
                    value = {repelem(obj.q_distribBin,obj.PulsesPerBin*obj.oversampling{1}),repelem(obj.q_distribBin,obj.PulsesPerBin*obj.oversampling{2})};
                end
            end
            % some error checks
           if ~iscell(value), error('q_distribFull must be a cell. Something went wrong.'); end
            refsz = [obj.TotalNoPulses*obj.oversampling{1},obj.TotalNoPulses*obj.oversampling{2}]; % reference size
            if ~isequal(cellfun(@(x) size(x,1),value),refsz)
                error('Error generating q_distribFull. Output must be [1,2] cell of cells [sz1,1] / [sz2,1] with sz1/2 = %s',num2str(refsz));
            end
        end
                
        %% SET METHODS
        
        function obj = set.bins(obj,value)
            % number of bins and type of binning (static number of pulses
            % or variable number of pulses per bin); static is legacy
            if ~isscalar(value{1})
                error('bins must be cell array {scalar,char}')
            end
            obj.bins = value; % write number of bins
            
            value{2} = lower(char(value{2}));
            allowed = {'static','variable'};
            if ismember(value{2},allowed)
                obj.bins{2} = value{2};
            else
                warning(strcat('bins{2} must be of type:',32,strjoin(allowed,', ')))
                fprintf('bins{2} defaulting to %s\n',allowed{1});
                obj.bins{2} = allowed{1};
            end
            obj = calcPulsesPerBin(obj); % update PulsesPerBin
        end
        
        function obj = set.precision(obj,value)
            % cell {calculation,storage}
            % calculation saves some speed doing FFTs, storage is worth
            % considering if lots of pulses are evaluated
            value = cellfun(@(x) lower(char(x)),value,'un',0);
            allowed = {'fp16','fp32','fp64'};
            if ismember(value{1},allowed(2:3))
                obj.precision(1) = value(1);
            else
                obj.precision(1) = allowed(3);
                warning(strcat('Precision for calculations must be of type:',32,strjoin(allowed(2:3),', '),', set to',32,obj.precision{1}))
            end
            if ismember(value{2},allowed)
                obj.precision(2) = value(2);
            else
                obj.precision(2) = allowed(2);
                warning(strcat('Precision for storage must be of type:',32,strjoin(allowed(1:3),', '),', set to',32,obj.precision{2}))
            end
        end
        
        function obj = set.Heatsource(obj,value)
            knownSources = {'point','gauss','ring'};
            if ismember(value,knownSources)
                obj.Heatsource = value;
            else
                error(strcat('Heatsource must of type:',32,strjoin(knownSources,', ')))
            end
        end
                
        function obj = set.EnergyRamp(obj,value)
            if ~isnumeric(value)
                error('EnergyRamp must be numeric')
            end
            if min(value(:)) < 0
                error('EnergyRamp must not have negative values')
            end
            if abs(1-max(value(:))) > eps
                obj.EnergyRamp = value./max(value(:));
                warning('Max value of provided EnergyRamp != 1 --> EnergyRamp has been normalized')
            else
                obj.EnergyRamp = value;
            end
        end
        
        function obj = set.oversampling(obj,value)
            obj.oversampling{1} = 1;
            % if value{2} is integer and >1 then we set it, else set to 1
            if (floor(value{2}) == value{2}) && value{2} >= 1
                obj.oversampling{2} = value{2};
            else
                obj.oversampling{2} = 1;
                warning('Oversampling set to 1. Oversampling must be integer between 1 and infinity')
            end
            obj = calcPulsesPerBin(obj); % update PulsesPerBin
        end
        
        function obj = set.TotalNoPulses(obj,value)
            if ~(2^nextpow2(value) == value)
                obj.TotalNoPulses = 2^nextpow2(value);
                warning(strcat('Total number of pulses is set to next power of 2:',32,num2str(obj.TotalNoPulses)))
            else
                obj.TotalNoPulses = value;
            end
            obj = calcPulsesPerBin(obj); % update PulsesPerBin
        end
        
        function obj = set.Depths_Eval(obj,value)
            % if a scalar or vector is provided write it to cell 1
            if isnumeric(value{1})
                obj.Depths_Eval{1} = value{1};
                obj.Depths_Eval{2} = 'userdefined';
            elseif isequal(size(value),[1,2])
                % it is a cell and must be of appropriate size
                obj.Depths_Eval = value;
            else
                error('Provided Depths_Eval does not match expectation')
            end
            % trigger generation of z eval positions
            obj = calcZevalPositions(obj);
        end
        
        function obj = set.Heatsource2Dsize(obj,value)
            % length of the Heatsource at any given pulse
            % cell containing 2D Matrix (TotalNoPulses,:) or array of cells (TotalNoPulses,:))
            % array size validation together with qdistrib verify_lheatqvect_cellarrays(obj)
            if isequal(size(value),[1,1]) && ~isscalar(value{:})
                value = num2cell(value{:},2); % if matrix then store as array of cells (TotalNoPulses,:))
            elseif size(value,2) > 1
                error('Heatsource2Dsize must be (1,1) cell containing 2D array, scalar, or (npulses,1) cell array of scalars / row vectors');
            end
            if ~any(cellfun(@isnumeric,value))
                error('Heatsource2Dsize must be cell array of numerics (scalars or vectors).')
            end
            for i = 1:numel(value)
                if any(value{i} < 0,'all')
                    error('Error writing Heatsource2Dsize; Size of heat source must not have negative values!')
                end
            end
            obj.Heatsource2Dsize = value; % write
        end        
        
        function obj = set.Heatsource2Dlen(obj,value)
            % length of the Heatsource at any given pulse
            % cell containing 2D Matrix (TotalNoPulses,:) or array of cells (TotalNoPulses,:))
            % array size validation together with qdistrib verify_lheatqvect_cellarrays(obj)
            if isequal(size(value),[1,1]) && isvector(value{:}) && ~isscalar(value{:})
                value = num2cell(value{:},2); % if matrix then store as array of cells (TotalNoPulses,:))
            end
            if ~any(cellfun(@isnumeric,value))
                error('Heatsource2Dlen must be cell array of numerics (scalars or vectors).')
            end
            for i = 1:numel(value)
                value{i}(value{i} == 0) = eps; % just to avoid singularities
                if any(value{i} < 0,'all')
                    error('Error writing Heatsource2Dlen; Length of heat source must not have negative values!')
                end
                if ~isscalar(value{i}) % then it must be a vector
                    diffval = diff(value{i}); % forward diff
                    if ~isempty(diffval) && any(diffval<0) % then vector is decreasing at some point
                        error('If Heatsource2Dlen vector defines heat source positions, input should be continually increasing!')
                    end
                end
            end
            obj.Heatsource2Dlen = value; % write
            obj = calcZevalPositions(obj); % now get z-Positions if requested
            obj = calcPulsesPerBin(obj); % update PulsesPerBin
        end
        
        function obj = set.CorrfactAbs(obj,value)
           % CorrfactAbs is meant to account for increased absorption, i.e.
           % through Gouffé or Raytracing. For each pulse (geometry) one
           % scalar should be provided, i.e. CorrfactAbs is a vector
           % (1,TotalNoPulses)
           if ~isvector(value)
               error('Energy correction factor CorrfactAbs must be a vector');
           end
           if min(value) < 1
               warning('CorrfactAbs contains values < 1. This means residual energies lower than pulse energy * eta_abs occur');
           end            
           obj.CorrfactAbs = value;
        end

        function obj = set.Q_distrib_2D(obj,value)
            % cell containing vector (1,:), string / char (1,1) or 2D Matrix (TotalNoPulses,:) or array of cells (TotalNoPulses,:)) 
            % for despription of energy distribution along 2D heat source from 3D point sources
            % note: don't handle input like set.Heatsource2Dlen bc we have different interpolation implementations for qdistribbin
            % (no preliminary conversion like num2cell(value{:},2)
            if any(cellfun(@isnumeric,value))
                % we can't have any zeros or in edge case we divide by zero
                for i = 1:numel(value)
                    value{i}(value{i} == 0) = eps;
                    if any(value{i} < 0,'all')
                        error('Error writing Q_distrib_2D; Energy distribution must not have negative values!')
                    end
                end
            end
            obj.Q_distrib_2D = value;
        end
    end
    
    % Private and Static methods
    methods (Access = private)
        function obj = calcZevalPositions(obj)
            % check if Depths_Eval is provided as cell(string,double)
            verify = cellfun(@class,obj.Depths_Eval,'un',0);
            if (isequal(verify,{'string','double'}) ||  isequal(verify,{'char','double'}))
                switch lower(obj.Depths_Eval{1})
                    % allowed syntax is {'determine',<value>} where value = number of evaluated positions
                    case 'determine'
                        if ~isempty(obj.Heatsource2Dlen)
                            % pick last depth from last element of last cell
                            % make sure we get one depth below the last
                            % heat input
                            % obj.z_eval_depths = linspace(0,obj.Heatsource2Dlen{end}(end),obj.Depths_Eval{2});
                            obj.z_eval_depths = linspace(0,max(cellfun(@max,obj.Heatsource2Dlen))*1.1,obj.Depths_Eval{2});
                        else % reset to zero
                            disp('Z-evaluation depths will be calculated once Heatsource2D lengths have been provided')
                            obj.z_eval_depths = 0;
                        end
                    otherwise
                        error(strcat('z-Position eval case',32,obj.Depths_Eval{1}),32,'not defined')
                end
            elseif isnumeric(obj.Depths_Eval{1}) && ~isempty(obj.Depths_Eval{1})
                % else it should be vector or scalar, use user provided values
                obj.z_eval_depths = obj.Depths_Eval{1};
            end
        end
        function obj = calcPulsesPerBin(obj)
            % Calculates property PulsesPerBin when Heatsource2Dlen
            % / TotalNoPulses / bins is written
            % PulsesPerBin implicitly or explicitly defines the pulses which are
            % evaluated with a shared heat kernel
            % (i.e. instead of one heat kernel per drill depth we use one heat kernel for PulsesPerBin)
            if ~any(cellfun(@isempty,obj.bins))
                if ismember(obj.bins{2},{'variable'}) && ~isempty(obj.Heatsource2Dlen) % either static or variable, requires Heatsource2Dlen
                    if max(cell2mat(cellfun(@numel,obj.Heatsource2Dlen,'un',0))) > 1
                        % default to static in case user provides explicit position and energy distribution for each pulse "Weber-edgecase"
                        warning('User-provided position vectors / energy distributions currently not implemented with variable bins. Defaulting to static binning.')
                        obj.bins{2} = 'static';
                    elseif (min(diff(cell2mat(obj.Heatsource2Dlen))) < 0) || (sum(diff(cell2mat(obj.Heatsource2Dlen))) == 0)
                        % else provided is cell array of scalars.
                        % if above statements are true then: either not steadily growing OR constant length
                        % in either case static binning is required. assumption is that user specifies temporally varying energy distribution
                        % hence binning cannot be derived from length of heat source and is instead applied after constant PulsesPerBin
                        % -- variable binning requires monotonously growing Heatsource2Dlen on linear grid --
                        obj.bins{2} = 'static';
                        warning('Heatsource2Dlen is not steadily growing. Cannot be drilling. Defaulting to static bining')
                    else
                        % assumption is that varying energy distributions correlate primarily with depth of borehole.
                        % distribute pulses per bin such that the
                        % heatkernels are evenly distributed opposed to every n pulses a new heatkernel. this should be
                        % beneficial because it allows for better resolution at the start of the drilling with
                        % binning, while at the end lots of pulses may share the same heatkernel when process stagnates
                        %
                        %              oversampling = 1            |             oversampling{2}
                        % ARRAY: idx_start idx_end mean_idx pulses | idx_start idx_end pulses oversampling_check
                        % mean_idx is used to extrapolate relevant unbinned values from binned
                        depths = cell2mat(obj.Heatsource2Dlen);
                        mindepth = min(depths); maxdepth = max(depths);
                        depthbins = linspace(mindepth,maxdepth,obj.bins{1});
                        val = ones(obj.bins{1},1); % holds indices where jump to next bin is required
                        for i = 2:length(depthbins)
                            val(i) = find(depths < depthbins(i),1,'last');
                            if val(i) < val(i-1); val(i) = val(i-1)+1; end
                        end
                        bin_idx = NaN(obj.bins{1},8); % init
                        bin_idx(:,1) = val; % startindex for each bin
                        bin_idx(:,2) = [val(2:end)-1; length(obj.Heatsource2Dlen)]; % endindex for each bin
                        bin_idx(:,3) = ceil((bin_idx(:,2)+bin_idx(:,1))/2); % mean index for each bin (use to construct relevant binned quantities)
                        bin_idx(:,4) = 1+bin_idx(:,2)-bin_idx(:,1); % number of pulses for that bin
                        bin_idx(:,6) = cumsum(bin_idx(:,4)*obj.oversampling{2}); % end indices oversampling{2}
                        bin_idx(:,5) = bin_idx(:,6)-(bin_idx(:,4)*obj.oversampling{2})+1; % start indices @ oversampling{2}
                        bin_idx(:,7) = 1+bin_idx(:,6)-bin_idx(:,5); % number of pulses / timesteps for that bin @ oversampling{2}
                        bin_idx(:,8) = bin_idx(:,7)./bin_idx(:,4); % error checking: must be NaN or oversampling{2}
                        
                        if (sum(bin_idx(:,4)) ~= obj.TotalNoPulses) || (sum(bin_idx(:,7)) ~= obj.TotalNoPulses*obj.oversampling{2})
                            error('calcPulsesPerBin (variable) resulted in incorrect total number of pulses / timesteps after binning')
                        end
                        obj.PulsesPerBin = bin_idx;
                    end
                else % default to static
                    obj.PulsesPerBin = obj.TotalNoPulses/obj.bins{1};
                end
            else
                obj.PulsesPerBin = NaN;
            end
        end
    end
    
    methods (Static, Access = private)
        function output = binqvectmatrix(input,bins,TotalNoPulses,PulsesPerBin)
            if strcmp(bins{2},'variable')
                output = input(PulsesPerBin(:,3),:); % just extract corresponding to mean idx
            else
                % generate binned array of q vectors for 3DLine energy distribution for convolution calc
                output = interp2(1:size(input,2), linspace(1,bins{1},TotalNoPulses), ... [Y old, X old, Input, Y new, X new]
                          input, 1:size(input,2), (1:bins{1}).');
            end
            output = num2cell(output,2); % convert to cell, we always want cell arrays of row vectors
        end
        function output = unbinqvectmatrix(input,bins,TotalNoPulses,PulsesPerBin,oversampling)
            output = {cell(TotalNoPulses*oversampling{1},1),cell(TotalNoPulses*oversampling{2},1)}; % init
            if strcmp(bins{2},'variable')
                for i = 1:bins{1}
                    if PulsesPerBin(i,4)
                        % then repeat the elements according to the indices stored in PulsesPerBin(i,4) which means
                        % indices obj.PulsesPerBin(i,1) through obj.PulsesPerBin(i,2) (obj.PulsesPerBin(i,4) pulses times oversampling 1)
                        output{1}(PulsesPerBin(i,1):PulsesPerBin(i,2),1) = num2cell(repmat(input(i,:),[PulsesPerBin(i,4)*oversampling{1},1]),2);
                        % indices obj.PulsesPerBin(i,5) through obj.PulsesPerBin(i,6) (obj.PulsesPerBin(i,4) pulses times oversampling 2)
                        output{2}(PulsesPerBin(i,5):PulsesPerBin(i,6),1) = num2cell(repmat(input(i,:),[PulsesPerBin(i,4)*oversampling{2},1]),2);
                    end
                end
            else
                % generate full arrays of q vectors for 3DLine energy distribution for explicit calc
                for i = 1:length(oversampling)
                    tmp = interp2(1:size(input,2), linspace(1,TotalNoPulses*oversampling{i},bins{1}), ... [Y old, X old, Input, Y new, X new]
                           input, 1:size(input,2), (1:TotalNoPulses*oversampling{i}).');
                    output{i} = num2cell(tmp,2); % convert to cell, we always want cell arrays of row vectors
                end
            end
        end
        function chirp_out = fix_reprates(reprate,chirp)
            % modifies an input chirp vector such that the closest
            % achievable reprate is calculated based on the maxiumum
            % provided reprate using integer dividers
            compvect = reprate./(1:100); % possible reprates, first 100 integer divider
            chirp_out = zeros(size(chirp));
            for i = 1:length(chirp)
                compare = repmat(chirp(i),[1 (length(compvect))])-compvect;
                [~,idx] = min(abs(compare));
                chirp_out(i) = compvect(idx);
            end
        end
        function output = generateID(NumOfDigits)
            output = char(java.util.UUID.randomUUID.toString);
            output = output(1:NumOfDigits);
            % gen random lowercase char 
            randstr = char(97:106); randidx = round(9*rand(1)+1,0);
            % append character in case we want to use this as variable name (first char digit not allowed)
            output = strcat(randstr(randidx),output);
        end
    end
    
    methods (Static)
        function useallcores()
            % 2020a forward seems to default to only logical cores
            % hower hyperthreading benefits the calculations so we enable
            % if possible
            try  %#ok<TRYNC>
                maxNumCompThreads(str2double(getenv('NUMBER_OF_PROCESSORS')));
            end
        end
        
        function obj = update_timesteps(obj,input,variant,reprate_laser)
            % <<< Timestep generation explicitly requires external input >>>
            %
            % Generate timestep and timestep_hq vector corresponding to laser
            % pulse to pulse time (i.e. frep^-1)
            % use obj = get_timesteps(obj,laser.frep); % where laser.frep
            % is current reprate and hence 1/laser.frep is pulse to pulse
            % length
            %
            % "input" must either be laser.frep (timestep / pulse to pulse length) for case 'constant' timestepping
            % or a user generated vector of length(1,TotalNoPulses) for case 'variable' for frequency chirps
            % chirps are currently not implemented for the heat accumulation calculation
            %
            % Deprecated: Generate timesteps directly
            % obj.timesteps = (laser_tpp : laser_tpp : obj.TotalNoPulses*laser_tpp);
            % obj.timesteps_hq = (laser_tpp/obj.oversampling{2} : laser_tpp/obj.oversampling{2} : obj.TotalNoPulses*laser_tpp);
            %
            % New Method: Generate index vector and oversampled index vector. If a
            % timestepping vector is provided an oversampled vector can be generated
            % using linear interpolarion
            
            if ismember(class(variant),{'char','string'})
                obj.timesteps{3} = char(lower(variant));
            else
                error('variant must be string or character "constant" or "variable"')
            end
            
            switch lower(variant)
                case 'constant'
                    % here we can generate timesteps directly and calculate
                    % timestepping from the constant reprate which is provided
                    if numel(input) ~= 1
                        error('laser frequency must be of size [1,1] and provided in Hz')
                    end
                    indexes = (1 : 1 : obj.TotalNoPulses);
                    indexes_hq = (1/obj.oversampling{2} : 1/obj.oversampling{2} : obj.TotalNoPulses);
                    obj.timesteps{1} = (1/input : 1/input : obj.TotalNoPulses*1/input);
                    obj.timesteps{2} = interp1(indexes,obj.timesteps{1},indexes_hq,'linear','extrap');
                    obj.reprates{1} = input*ones(1,obj.TotalNoPulses);
                    obj.timestepping{1} = 1./obj.reprates{1};
                    obj.timestepping{2} = repelem(obj.timestepping{1},obj.oversampling{2})./obj.oversampling{2};
                    obj.reprate_max = input;
                case {'variable','variable-intdivider'}
                    % here the input is a frequency chirp vector for the
                    % laser settings: 1/chirp yields timings between pulses
                    % and its cumulative sum are the evaluation times (the
                    % timesteps variable). Calculation order is reversed in
                    % this case.
                    %
                    % chirp is more involved, each heat input requires its
                    % own temporally shifted time vector, i.e. its solution
                    % T(t,r) is placed into the superposition matrix at the
                    % index when its heat input occurs
                    % this must be evaluated explicity, convolution method
                    % is not applicable
                    if length(input) == obj.TotalNoPulses
                        if strcmpi(variant,'variable-intdivider')
                            if ~exist('reprate_laser','var')
                                error('timesteps "variable-intdivider": Maximum laser pulse repetition rate as last argument required.')
                            end
                            % only use reprates with integer dividers
                            obj.reprates{1} = obj.fix_reprates(reprate_laser,input);
                        else
                            obj.reprates{1} = input;
                        end
                        obj.timestepping{1} = 1./obj.reprates{1};
                        obj.timestepping{2} = repelem(obj.timestepping{1},obj.oversampling{2})./obj.oversampling{2};
                    else
                        error('frequency chirp must be size [1,TotalNoPulses]');
                    end
                    obj.timesteps{1} = cumsum(obj.timestepping{1});
                    obj.timesteps{2} = cumsum(obj.timestepping{2});
                    % if timestepping is variable we can't use convolution
                    % method anyway and therefore there is no benefit to
                    % binning (need to calculate for each pulse anyway)
                    if obj.bins{1} < obj.TotalNoPulses
                        warning('variable timestepping: Binning has no computational benefit -> disabled')
                        obj.bins = {obj.TotalNoPulses,'static'};
                    end
                    obj.reprate_max = max(input);
                otherwise
                    error('timestep variant unknown [known: "variable" or "constant"')
            end
            obj.reprates{2} = repelem(obj.reprates{1},obj.oversampling{2});
        end
        function verify_lheatqvect_cellarrays(obj)
            if ~isempty(obj.Heatsource2Dlen) && ~isempty(obj.Q_distrib_2D)
                if max(cellfun(@(x) size(x,2), obj.Heatsource2Dlen)) > 1
                    % then Heatsource2Dlen is array of vectors and thus they
                    % describe the actual positions of the heat sources
                    % in this case corresponding vectors in Heatsource2Dlen and Q_distrib_2D
                    % must have equal number of elements or there is an error
                    verification = cellfun(@isequal, ...
                        cellfun(@size,obj.Heatsource2Dlen,'UniformOutput',false),...
                        cellfun(@size,obj.Q_distrib_2D,'UniformOutput',false));
                    bad_indexes = find(~verification).';
                    if ~isempty(bad_indexes)
                        error('Vectors for positions of heat sources and heat quantity do not match for indices:\n%s',...
                            mat2str(bad_indexes))
                    end
                end
            else
                warning('call to verify_lheatqvect_cellarrays while either Heatsource2Dlen or Q_distrib_2D is empty')
            end
            if obj.debug, fprintf('Consistency check of Heatsource2Dlen and Q_distrib_2D --- OK\n'); end
        end
    end
end
