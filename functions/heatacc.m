function results = heatacc(laser,mat,sim)
% This function calculates heat accumulation according to superposition of point sources for the half
% infinite body (including mirror sources about z = 0)
% The length of the heat source is calculated from the drill depth progress calculation where
% l_heat = l_drill(n_pulse)
% In contrast to the fundamental 2D solution there is no singularity 1/l_heat for
% l_heat -> 0. In addition, the relative energy distribution along the
% point sources distributed along the line is variable (q_vect in class optsheatsource)

%% Select appropriate settings if we use oversampling

% Relevant depdendence on oversampling here: timesteps, timestepping, HeatInputsPerBin
if sim.oversampling{2} ~= 1
    idx = [1,2];
else
    idx = 1;
end

%% Do we have constant or variable (chirped) timestepping?
% constant -> convolution, variable -> explicit
% debug_override -> both and compare
variant = sim.timesteps{3};
% variant = 'debug_override';

%% Information / Warnings
if strcmpi(sim.timesteps{3},'variable')
    warning('frequency chirped heat accumulation is significantly slower to calculate')
end

%% calculate drill depth in microns and absorption correction factor according to Gouffé
if isempty(sim.Heatsource2Dlen)
    disp('Calculating drilling progress...')
        [depth_array,Abs_Fact,Gouffe_Abs,drilldepth_max,energy_distribution,ablation_crater_radius] = drilldepth(laser,mat,sim); % get drilldepht and absorption correction factor
    if ~any(isnan(energy_distribution(:)))
        % sim.Q_distrib_2D = {energy_distribution}; % alternative pass as 1x1 cell containing matrix [TotalNoPulses,whatever]
        sim.Q_distrib_2D = num2cell(energy_distribution,2); % alternative pass as cell array
    end
    % write Heatsource2D lengths (for analytical 2D this should be corrected by lheat thermal
    sim.Heatsource2Dlen = {depth_array}; % pass as vector in cell or cell for each pulse, i.e. = num2cell(depth_array,2);
    sim.CorrfactAbs = Abs_Fact;
    sim.isdrilling = [1,ablation_crater_radius]; % set drill flag
else
    disp('Heatsource length as function of pulses has been provided...')
    % dirty hacks so I don't need case handling for standard plots
    depth_array = cell2mat(sim.Heatsource2Dlen); % 16.04.21 deprecated? shouldnt this be cell?
    Gouffe_Abs = ones(1,length(depth_array));
    drilldepth_max = max(depth_array(:));
end

% check if Q_distrib_2D and Heatsource2Dlen match (only applies to the case
% where Heatsource2Dlen is provided as vectors and explicitly defines positions of the heat sources
sim.verify_lheatqvect_cellarrays(sim)

%% store computationally expensive dependent variables locally
% store l_heat and q_vect for convolution based variant
loc_heatsource2DlenBin = sim.heatsource2DlenBin;
loc_heatsource2DsizeBin = sim.heatsource2DsizeBin;
loc_q_distribBin = sim.q_distribBin;

% store l_heat and q_vect for explicit calculation based variant
loc_heatsource2DlenFull = sim.heatsource2DlenFull;
loc_heatsource2DsizeFull= sim.heatsource2DsizeFull;
loc_q_distribFull = sim.q_distribFull;
loc_corrfactFull = sim.corrfactFull;

%% preallocate / initialize result arrays
numericalHeat_vertsum_z = initresultarray(sim);

% initialize storage for output q_vect (energy distribution along l_heat)
q_vect = cell(sim.bins{1},1); % output from heatkernel calc, atm only used for console output of totalnumsources

%% Make all following calculations normal and oversampled if requested
ttotal = tic;
for i = idx
    % constant timestepping (fixed reprate)
    if ismember(variant,{'constant','debug_override'})
        disp('Constant timestepping: Linear convolution through FFTs')
        % construct heat input / pulse energy matrix
        pulsematrix = genpulsematrix(sim,laser,mat,i,'new'); % 'legacy' for old kron variant        
        % Determine required zero-padding circconv -> linconv + pad to nextpow2 [a+b-1]
        zeropad = 2^nextpow2(length(sim.timesteps{i})+size(pulsematrix,2)-1);
        % determine limiting size of array dimension as nextpow2 to keep within memory boundaries [1Gb array size limit]
        len = getArraySize(length(sim.timesteps{i}),sim.bins{1},sim.precision{1});
        % determine indexing for overlap-add loop (if len < sim.bins{1} then process in chunks)
        lenvect = [0 len:len:sim.bins{1}];
        % start timer
        tstart = tic; 
        % count finished parameters
        counter = 1;
        % loop over all r positions
        for m = 1:length(sim.radial_position)
            sim.HeatsourceSettings.x_pos = sim.radial_position(m);
            % loop over all z positions
            for n = 1:length(sim.z_eval_depths)
                % select current z position
                sim.HeatsourceSettings.z_pos = sim.z_eval_depths(n);
                % reset / initialize result vector numericalHeat
                numericalHeat = zeros(1,length(sim.timesteps{i}));
                for k = 1:(sim.bins{1}/len)
                    % current overlap-add range
                    range = lenvect(k)+1:lenvect(k+1);
                    % generate heatkernels for current chunk
                    sim.HeatsourceSettings.l_heat = loc_heatsource2DlenBin(range);
                    % set current scalar size for heat source (scalar) or vector
                    sim.HeatsourceSettings.sz_heat = loc_heatsource2DsizeBin(range);
                    % set current q_vect
                    sim.HeatsourceSettings.q_vect = loc_q_distribBin(range);
                    % calculate all kernels for the specified heat source lenghts
                    [heatkernel_3DLine_array, q_vect_output] = get2D_Distributed_Heatkernels(sim.HeatsourceSettings,mat,sim.timesteps{i});
                    % convolve pulsematrix with heatkernels
                    switch sim.precision{1}
                        case 'fp32' % note: In testing always faster than dircrete conv, difference to fp64 negligible
                            overlapaddHeat = ifft(fft(single(full(pulsematrix(range,:))),zeropad,2).* ...
                                                  fft(single(heatkernel_3DLine_array),zeropad,2),[],2,'symmetric');
                        case 'fp64'
                            overlapaddHeat = ifft(fft(full(pulsematrix(range,:)),zeropad,2).* ...
                                                  fft(heatkernel_3DLine_array,zeropad,2),[],2,'symmetric');
                    end
                    % remove residual fft artifacts
                    overlapaddHeat(overlapaddHeat < 1e-12) = 0;
                    % remove values outside of kernel overlap
                    overlapaddHeat = overlapaddHeat(:,1:length(sim.timesteps{i})-1);
                    % superimpose / sum orthogonal to temporal dimension
                    numericalHeat(2:end) = numericalHeat(2:end)+sum(overlapaddHeat,1);
                    % console output
                    if (length(lenvect) > 2), fprintf('chunk %i/%i\n',k,length(lenvect)-1); end
                    if (n == 1 && i == 1), q_vect(range) = q_vect_output; end
                end
                % save result to cell array
                numericalHeat_vertsum_z{0+i,m}(n,:) = numericalHeat;
                % console output
                if len < sim.bins{1} || length(sim.z_eval_depths) > 1
                    current_time = toc(tstart); time_remain = round((current_time/counter)*(length(sim.radial_position)*length(sim.z_eval_depths))-current_time,0);
                    counter = counter+1; % advance
                    fprintf('est. remain: %3.2f [min]\n',time_remain/60)
                end
            end
        end
    end
    
	% variable timestepping - for frequency chirps (variable reprates)
    if ismember(variant,{'variable','variable-intdivider','debug_override'})
        disp('Variable timestepping: Explicit-method')
        % we need pulsevect for variable timestepping / explicit calculation
        pulsevect = sim.HeatInputs{i}.*sim.EnergyRampFull{i};
        % start timer
        tstart = tic;
        % count finished parameters
        counter = 1;
        % loop over all r positions
        for m = 1:length(sim.radial_position)
            sim.HeatsourceSettings.x_pos = sim.radial_position(m);
            % loop over all z positions
            for n = 1:length(sim.z_eval_depths)
                fprintf('[z,r] = %4.1f / %3.1f µm: ',sim.z_eval_depths(n)*1e6,sim.radial_position(m)*1e6)
                % select current z position
                sim.HeatsourceSettings.z_pos = sim.z_eval_depths(n);
                % reset / initialize result vector numericalHeat
                numericalHeat = zeros(1,sim.TotalNoPulses*sim.oversampling{i});
                for k = 1:sim.TotalNoPulses*sim.oversampling{i}
                    % if pulsevect(i) is zero we can skip
                    if pulsevect(k) ~= 0
                        % get timesteps for the current pulse [inf(1,index-1) cumsum(timestepping(index:end))];
                        current_timesteps = gettimesteps(k,sim.timestepping{i});
                        % set current Heat source length
                        sim.HeatsourceSettings.l_heat = loc_heatsource2DlenFull{i}(k);
                        % set current scalar beam radius for heat source (scalar) or vector
                        sim.HeatsourceSettings.sz_heat = loc_heatsource2DsizeFull{i}(k);
                        % set current q vect
                        sim.HeatsourceSettings.q_vect = loc_q_distribFull{i}(k);
                        % calculate kernels for the specified heat source length
                        [heatkernel_3DLine, q_vect_output] = get2D_Distributed_Heatkernels(sim.HeatsourceSettings,mat,current_timesteps);
                        % and scale with energy input times correction factor
                        tempfield_pulsek = (pulsevect(k)*laser.Ep*mat.eta_abs).* ...
                                           (loc_corrfactFull{i}(k)).*heatkernel_3DLine;
                        % optimization note: multiplying with infs is computationally efficient
                        % superimpose temp fields
                        numericalHeat = numericalHeat+tempfield_pulsek;
                        if (n == 1 && i == 1), q_vect(k) = q_vect_output; end
                    end
                    % console output
                    if mod(k,floor(length(pulsevect)/10)) == 0
                        fprintf('%1.0f%% - ',100*k/length(pulsevect));
                    end
                end
                % timeshift and save to cell array
                numericalHeat_vertsum_z{2+i,m}(n,:) = [0 numericalHeat(1:end-1)];
                % console output
                if length(sim.z_eval_depths) > 1
                    current_time = toc(tstart); time_remain = round((current_time/counter)*(length(sim.radial_position)*length(sim.z_eval_depths))-current_time,0);
                    counter = counter+1; % advance
                    fprintf('est. remain: %3.2f [min]\n',time_remain/60)
                end
            end
        end
    end
end
totalsimtime = toc(ttotal);
fprintf('Total number of superimposed unique heat sources (excluding mirror sources): %1.2e\n',sum(cellfun(@length,q_vect)))

%% Output complete results
results = drillheataccResults(); % init
results.laser = laser; results.sim = sim; results.mat = mat;
results.totalsimtime = totalsimtime;
results.qvect_kernel = q_vect;
results.drilldepth_max = drilldepth_max; % förster limit, req. for plot
results.gouffe = Gouffe_Abs; % required 4 plot 1

% dependent on calculation scheme some of the cells are empty, remove these write to results
if ismember(variant,{'constant'})
    results.TempResults_z = numericalHeat_vertsum_z(1:2,:); % saving convolution 
elseif ismember(variant,{'variable','variable-intdivider'})
    results.TempResults_z = numericalHeat_vertsum_z(3:4,:); % saving explicit
else
    error('Debug-Override is active. Output of results not defined / intended')
end

%% Plot in case of debug override
switch variant
    case 'debug_override'
        warning(strcat('Debug override active in',32,mfilename))
end

% figure; plot(numericalHeat_vertsum{1,1}); hold on, plot(numericalHeat_vertsum{2,1},'--k');
% figure; plot(numericalHeat_vertsum{1,1}-numericalHeat_vertsum{2,1})
% title('Difference between convolution and explicit result')

end

function output = gettimesteps(index,timestepping)
% outputs the correct evaluation timesteps for the pulse applied at [index]
% given the timestepping of the process
% setting to inf before the pulse exists ensures T = 0 at these times
% tested against anon and direct evaluation...this is fastest
% after profiling localfun > nofun > inline/anonymous
% timesteps_ForPulseN = @(index,timestepping) [inf(1,index-1) cumsum(timestepping(index:end))];
%
% optimization potential is given by computing a "timestepping" matrix and
% then getting all the heatkernels for the pulses at once...tbd
%
output = [inf(1,index-1) cumsum(timestepping(index:end))];
end

function output = initresultarray(sim)
% this function determines the required output cell array and preallocates
% the arrays within the cell array according to the precision specification
% this is done with the intent to avoid RAM out of memory during the computation
% because of excessive input parameters (otherwise this would occur slowly
% as the results are stored)
% further this ensures that the RAM-check for the convolutions is accurate
%
% if the host is a Windows PC a warning will be prompted for arrays that
% exceed 80% of the available RAM
% if the array requires more than 80% precision will be lowered from
% fp64 -> fp32 or fp32 -> fp16 to save space.
% MAC/Unix: A warning will be prompted for resulting arrays > 5Gb

% init output cell array
skip = false;
output = cell(4,length(sim.radial_position));
% determine sizes of output arrays in the cell array
arraysz = [length(sim.z_eval_depths),sim.TotalNoPulses*sim.oversampling{1},...
           length(sim.z_eval_depths),sim.TotalNoPulses*sim.oversampling{2}];
numofelems = [arraysz(1)*arraysz(2)*length(sim.radial_position),...
              arraysz(3)*arraysz(4)*length(sim.radial_position)];
% determine where to preallocate (debug override ignored)
if ismember(sim.timesteps{3},{'constant'})
    if sim.oversampling{2} == 1; indices = 1; % fill cell(1,:)
    else; indices = [1,2]; % fill cell(1:2,:)
    end
elseif ismember(sim.timesteps{3},{'variable','variable-intdivider'})
    if sim.oversampling{2} == 1; indices = 3; % fill cell(3,:)
    else; indices = [3,4]; % fill cell(3:4,:)
    end
else
    skip = true; % debug enabled, just return unfilled cell array
end
numofelems = numofelems(1:length(indices)); % sim.oversampling{2} == 1 then drop numofelems(2)
if ~skip
    if ispc % windows host
        sizelimit = memory; % get system RAM info
        %   in case of fast convolution array sizes are *2*2 (2 x fp64 = fp128, 2x signal and kernel, 2x zpadding)
        sizelimit = 0.8*sizelimit.MaxPossibleArrayBytes; % 80% of adressable RAM for results as limit
    else % OSX or Unix host
        sizelimit = 5e9; % 5gb
    end
    usr_precision = str2double(sim.precision{2}(3:4)); % extract "32" from 'fp32' etc
    allowed_precision = [2^6,2^5,2^4]; % fp64-fp32-fp16 | double-single-half
    % all sizes that might result dependent on precision
    size_possible = 2*sum(numofelems)*allowed_precision/8; % convert to byte (*1/8), overhead ~2x (factor 2)
    % check if any precision allows for storage within size constraint
    idx = find((sizelimit-size_possible)>0,1,'first');
    if ~isempty(idx)
        if allowed_precision(idx) >= usr_precision
            % user specified precision fits within constraint, don't comment
            output_precision = usr_precision;
            outputsz = [size_possible(allowed_precision == usr_precision)*1e-9,size_possible(allowed_precision == usr_precision)*1e-6];
            fprintf('Size of results will be > %2.3f Gb | %2.3f Mb.\n',outputsz(1),outputsz(2))
        else
            output_precision = allowed_precision(idx);
            warning('Numeric precision of output has been set to %i bit in order to avoid RAM Out-Of-Memory.',output_precision)
            warning('Size of results will be > %2.2f Gb | %2.3f Mb.\n',size_possible(idx)*1e-9,size_possible(idx)*1e-6)
        end
    elseif ~ispc
        output_precision = allowed_precision(3);
        warning('Results exceed at least 5Gb. Cannot verify on Mac/UNIX. Check RAM usage. Results precision set to %i bit.',output_precision)
        warning('Size of results will be > %2.2f Gb | %2.3f Mb.\n',size_possible(3)*1e-9,size_possible(3)*1e-6)
    else
        error('Size of results exceed %2.2f Gb EVEN at 16 bit precision. Reduce complexity / paramaters / get more RAM.\n',size_possible(3)*1e-9)
    end
    % populate cells with specified precision
    for i = 1:length(sim.radial_position)
        if output_precision == 64
            output{indices(1),i} = double(zeros(arraysz(1),arraysz(2)));
            if length(indices) == 2
                output{indices(2),i} = double(zeros(arraysz(3),arraysz(4)));
            end
        elseif output_precision == 32
            output{indices(1),i} = single(zeros(arraysz(1),arraysz(2)));
            if length(indices) == 2
                output{indices(2),i} = single(zeros(arraysz(3),arraysz(4)));
            end
        elseif output_precision == 16
            output{indices(1),i} = half(zeros(arraysz(1),arraysz(2)));
            if length(indices) == 2
                output{indices(2),i} = half(zeros(arraysz(3),arraysz(4)));
            end
        end
    end
end
end

function pulsematrix_out = genpulsematrix(sim,laser,mat,index,implementation)
% pulsematrix is a sparse block diagonal 2d matrix which contains the residual
% energies that apply to their respective binned heatkernel
% the following implementations are nearly identical, however the variant using
% kron cannot resolve variations WITHIN a bin as accurate, i.e. pulse energy ramps or
% on/off cycles may either get resolved badly or not at all using the
% kronecker product variant; kronecker is LEGACY implementation
implementation = lower(char(implementation)); % parse implementation (legacy forces kron if possible)
if nargin < 5 % assume new
    implementation = 'new';
end
if strcmp(implementation,'legacy') && strcmp(sim.bins{2},'variable')
   warning('Legacy (kron) pulsematrix not compatible with variable bins. Defaulting to generalized sparse block diagonal matrix.')
   implementation = 'new';
end
switch implementation
    case 'legacy'
        % legacy; ONLY WITH BINS 'static' & 'legacy'
        pulsevect_src = sim.EnergyRampBin.*laser.Ep*mat.eta_abs.*sim.corrfactBin;
        if ~isequal(size(pulsevect_src),[1,sim.bins{1}])
            error('pulsevect is not 1D row vector. Check input.');
        end
        pulsematrix_out = kron(spdiag(pulsevect_src),sim.HeatInputsPerBin{index});
    otherwise
        % sparse blockdiag - conserves energy for arbitrary binning and can handle arbitrary on/off cycles
        pulsevect_src = laser.Ep*mat.eta_abs.*sim.HeatInputs{index}.*sim.EnergyRampFull{index}.*sim.corrfactFull{index};
        if ~isequal(size(pulsevect_src),[1,sim.TotalNoPulses])
            error('pulsevect is not 1D row vector. Check input.');
        end
        if strcmp(sim.bins{2},'variable')
            pulsematrix = cell(sim.bins{1},1); rows = 1; % init
            for i = 1:sim.bins{1}
                if sim.PulsesPerBin(i,4)
                    % then repeat the elements according to the indices stored in PulsesPerBin(i,4) which means
                    if index == 1 % oversampling{1}
                        % indices obj.PulsesPerBin(i,1) through obj.PulsesPerBin(i,2) (obj.PulsesPerBin(i,4) pulses times oversampling 1)
                        len = length(pulsevect_src(sim.PulsesPerBin(i,1):sim.PulsesPerBin(i,2)));
                        pulsematrix{i} = sparse(zeros(rows,len));
                        pulsematrix{i}(rows,:) = sparse(pulsevect_src(sim.PulsesPerBin(i,1):sim.PulsesPerBin(i,2)));
                        rows = 1; % reset
                    elseif index == 2 % oversampling{2}
                        % indices obj.PulsesPerBin(i,5) through obj.PulsesPerBin(i,6) (obj.PulsesPerBin(i,4) pulses times oversampling 2)
                        len = length(pulsevect_src(sim.PulsesPerBin(i,5):sim.PulsesPerBin(i,6)));
                        pulsematrix{i} = sparse(zeros(rows,len));
                        pulsematrix{i}(rows,:) = sparse(pulsevect_src(sim.PulsesPerBin(i,5):sim.PulsesPerBin(i,6)));
                        rows = 1; % reset
                    else
                        error('This case should not exist.')
                    end
                else
                    rows = rows+1; % add row of zeros if no pulse is applied in this bin
                end
            end
        else
            % else sparse blkdiag static binning
            pulsematrix = reshape(pulsevect_src,[],sim.bins{1}).';
            pulsematrix = num2cell(sparse(pulsematrix),2);
        end
        pulsematrix_out = blkdiag(pulsematrix{:}); % note that empty cells are ignored (hence the 'rows' handling)
end
% verify array size is correct
arraysz = [size(pulsematrix_out);sim.bins{1},sim.TotalNoPulses*sim.oversampling{index}]; % [is, is; ref, ref]
if ~isequal(arraysz(1,:),arraysz(2,:))
    error('genpulsematrix: Pulsematrix is of size %ix%i but it should be of size %ix%i\n',...
        arraysz(1,1),arraysz(1,2),arraysz(2,1),arraysz(2,2))
end
% output energy, highlights conceptual problems with legacy kron and heavy binning
% conservation of energy should for arbitrary binning with spblkdiag
switch implementation
    case 'legacy'
        fprintf('Total energy static kron pulsematrix: %3.5f Joule --- ',full(sum(pulsematrix_out(:))))
    otherwise
        if strcmp(sim.bins{2},'variable')
            fprintf('Total energy variable blkdiag pulsematrix: %3.5f Joule --- ',full(sum(pulsematrix_out(:))))
        else
            fprintf('Total energy static blkdiag pulsematrix: %3.5f Joule --- ',full(sum(pulsematrix_out(:))))
        end
end
Ein_total = sum(laser.Ep*sim.HeatInputs{1}.*sim.EnergyRampFull{1});
fprintf('Total energy in: %3.5f Joule\n',Ein_total);
% user could provide inputs such that residual energy exceeds maximum pulse energy. throw error
peakresEnergy = 100*max(pulsevect_src)./(laser.Ep); % maximum applied residual energy [%], cannot exceed nomimal max laser pulse energy
fprintf('Peak residual energy is %2.2f%% of nominal maximum irradiated pulse energy\n',peakresEnergy);
if peakresEnergy > 100 %#ok<BDSCI>
    error('Peak residual energy cannot exceed nomimal maximum irradiated pulse energy. Bad user input.');
end
end