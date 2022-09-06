classdef drillheataccResults
% Containts results from drill heat accumulation calculation and calculates
% some useful metrics for plotting and analysis
    properties
        UUID                     char
        totalsimtime             double = NaN            
        EnableInterp             logical = false % this is only for interpolation of q_vect and heatsource2D length
                                                 % it serves only a benefit for visualization and may consume quite a bit of time
        interpVal(1,1)           double = 20 % for plotting, applies to q_vect_interp and heat2Dlen_interp
        laser                    initlaser
        mat                      materialprops
        sim                      simopts
        qvect_kernel             cell
        drilldepth_max           double
        gouffe                   double
        TempResults_z            cell
        MaxTemp_z                double
    end
    properties (SetAccess = protected)
        q_vect_interp            double
        heat2Dlen_interp         double
    end
    properties (Hidden, SetAccess = protected)
        interp_done(1,2)         cell = {false,NaN} % (logical / executed interpval)
    end
    properties (Dependent)
        TempResults_zrt             double % reshapes results into array [depth,radius,time/pulses];
        AvgLaserPower               double
        E_res                       double % effective residual energy in Joule
        E_in                        double
        MeanmaxavgLaserPower(1,2)   cell
        Qtotal_in(1,1)              double
        MillimeterPerJoule(1,1)     double
    end
    
    methods
                
        % Constructor
        function obj = drillheataccResults()
            %UNTITLED Construct an instance of this class
        end
        
        % Dependent
        function value = get.UUID(obj)
            % Calculate averager laser power Vector
            value = obj.sim.UUID;
        end
                
        function value = get.TempResults_zrt(obj)
            % reshapes results into array [depth,radius,time/pulses] upon request
            arraysz = [obj.sim.TotalNoPulses,length(obj.sim.z_eval_depths),length(obj.sim.radial_position)];
            % init zrt array
            switch obj.sim.precision{2}
                case 'fp64'
                    value = double(nan(circshift(arraysz,2)));
                case 'fp32'
                    value = single(nan(circshift(arraysz,2)));
                case 'fp16'
                    value = half(nan(circshift(arraysz,2)));
                otherwise % default to single
                    value = single(nan(circshift(arraysz,2)));
            end
            for i = 1:arraysz(3)
                value(:,i,:) = reshape(obj.TempResults_z{1,i},[arraysz(2),1,arraysz(1)]);
            end
        end
        
        function value = get.AvgLaserPower(obj)
            % Calculate averager laser power Vector
            value = (obj.sim.EnergyRampFull{1}.*obj.sim.HeatInputs{1}.*obj.laser.Ep)./ ...
                obj.sim.timestepping{1};
        end
        
        function value = get.MaxTemp_z(obj)
            % collapse temporal dimension -> max temp at a given depth/radius which occurs during process
            szz = length(obj.sim.z_eval_depths);
            szr = length(obj.sim.radial_position);
            value = NaN(szz,szr);
            for i = 1:szr
                for j = 1:szz
                    value(j,i) = max(obj.TempResults_z{1,i}(j,:));
                end
            end
        end
        
        function value = get.MeanmaxavgLaserPower(obj)
            % total energy incident divided by complete time until time(end of process)
            value{1} = sum(obj.sim.EnergyRampFull{1}.*obj.sim.HeatInputs{1}.*obj.laser.Ep)/obj.sim.timesteps{1}(obj.sim.PulseRange(2));
            value{2} = max(obj.AvgLaserPower);
        end
        
        function value = get.Qtotal_in(obj)
            % sum of total energy incident
            value = sum(obj.sim.EnergyRampFull{1}.*obj.sim.HeatInputs{1}.*obj.laser.Ep);
        end
        
        function value = get.E_res(obj)
            value = obj.sim.EnergyRampFull{1}.*obj.sim.HeatInputs{1}.*obj.laser.Ep.*obj.sim.corrfactFull{1}.*obj.mat.eta_abs;
        end
        
        function value = get.E_in(obj)
            value = obj.sim.EnergyRampFull{1}.*obj.sim.HeatInputs{1}.*obj.laser.Ep;
        end
        
        function value = get.MillimeterPerJoule(obj)
            maxdepth = max(cellfun(@max,obj.sim.Heatsource2Dlen));
            value = (maxdepth*1e3)/obj.Qtotal_in;
        end
                
        function obj = set.interpVal(obj,value)
            obj.interpVal = value;
            disp('For execution remember to set EnableInterp to true')
            obj = CalcInterpValues(obj);
        end
        
        function obj = set.EnableInterp(obj,value)
            obj.EnableInterp = value;
            obj = CalcInterpValues(obj);
        end
                
        function obj = CalcInterpValues(obj)
            if obj.EnableInterp == true
                if obj.interp_done{1} == false || (obj.interp_done{1} == true && obj.interp_done{2} ~= obj.interpVal)
                    tinterp = tic; disp('Qvect Interpolation...')
                    len = (obj.sim.TotalNoPulses*obj.sim.oversampling{1}); % array len / timesteps
                    if ~strcmpi('constant',obj.sim.Q_distrib_2D)
                    % interpolate qdistribbin to vector of length interpval
                    qdistribin = cellfun(@(x) interp1(1:length(x),x,linspace(1,length(x),obj.interpVal)),obj.sim.q_distribBin,'un',0);
                    % normalize
                    qdistribin = cellfun(@(x) x./sum(x), qdistribin,'un',0);
                    % remove nans if exist
                    qdistribin = obj.fixnan(qdistribin);
                    % generate qdistribfull for plotting from interpolated binned variant
                    if strcmp(obj.sim.bins{2},'variable')
                        value = cell(len,1);
                        for i = 1:obj.sim.bins{1}
                            if obj.sim.PulsesPerBin(i,4)
                                % then repeat the elements according to the indices stored in PulsesPerBin(i,4) which means
                                % indices obj.PulsesPerBin(i,1) through obj.PulsesPerBin(i,2) (obj.PulsesPerBin(i,4) pulses times oversampling 1)
                                value(obj.sim.PulsesPerBin(i,1):obj.sim.PulsesPerBin(i,2),1) = repelem(qdistribin(i),obj.sim.PulsesPerBin(i,4)*obj.sim.oversampling{1}).';
                            end
                        end
                    else
                        value = repelem(qdistribin,obj.sim.PulsesPerBin*obj.sim.oversampling{1});
                    end
                    else
                        % if energy distrib = constant just fill with ones where sum = 1
                        value = num2cell(repmat(ones(1,obj.interpVal)./obj.interpVal,[len,1]),2);
                    end
                    % q_vect does not know if energy = 0; fix these by setting to zeros of length interpval
                    energyrampfull_loc = obj.sim.EnergyRampFull{1}; % avoid repeated calls to dependent var
                    for i = 1:(obj.sim.TotalNoPulses*obj.sim.oversampling{1})
                        if energyrampfull_loc(i) == 0
                            value(i) = {NaN(1,obj.interpVal)}; % interpolation later requires minimum of 2 points
                        end
                    end
                    obj.q_vect_interp = cell2mat(value); % convert to 2d array
                    % now generate x-axis positions for the heat sources
                    disp('Heatsource2Dlen Interpolation...')
                    obj.heat2Dlen_interp = zeros(len,obj.interpVal); % init
                    depths_loc = obj.sim.Heatsource2Dlen; % drill depths for / after each pulse
                    % atm warn that qvect interp is not implemented for user provided explicit q position / hs len (only linear grid currently)
                    % implementation requires if depths_loc(i) = vector then interpolate nonlinear qvect/hslen grid to linear grids
                    warnenable = true;
                    for i = 1:len
                        if warnenable && ~isscalar(depths_loc{i})
                            warning(sprintf(['Visualization of energy distribution is only implemented for linear grid\n',...
                                'This warning is prompted because user-defined nonlinear energy distribution\n',...
                                'and / or user-defined heatsource positions were detected.\n',...
                                'Calculation is correct, but energy "dots" visualization is not.'])); %#ok<SPWRN>
                            warnenable = false; % disable further outputs
                        end
                        obj.heat2Dlen_interp(i,:) = linspace(0,depths_loc{i}(end),obj.interpVal);
                    end
                    % disable Interpolation flag and write
                    obj.EnableInterp = false;
                    obj.interp_done = {true,obj.interpVal};
                    fprintf('Interpolation computed after %3.2fs\n',toc(tinterp));
                end
            end
        end
        
        %% plots defined in external files for easier debugging
        plot1(obj) % metrics
        plot2(obj) % fixed z,r over t
        plot3(obj) % fixed r,t over z
        plotcfg = plot4(obj,value) % all r,z over t; 2nd argument is plotcfg which may be passed
        % -> simply call results.plot4(plotcfg);
        
        function plotcfg = plotall(obj)
            % simply calls all plots
            % e.g. call: plotcfg = results.plotall;
            obj.plot3
            obj.plot2
            obj.plot1
            plotcfg = obj.plot4;
        end
    end
    
    methods (Static)
        function x = fixnan(x)
            % recursively replace zeros by eps for numerics / cells
            if isnumeric(x)
                x(isnan(x)) = eps;
            elseif iscell(x)
                x = cellfun(@drillheataccResults.fixnan,x,'un',0);
            end
        end
    end
        
    methods (Static, Access = private)
        % note plot4gui / plot4update are defined as methods here bc they
        % need to scope each other, i.e. call each other via
        % classname.method
        output = plot4gui(referencefig,uuid,lastframe,initylims);
        [] = plot4update(data);
        % / plot4gui, plot4update
    end
    
end

