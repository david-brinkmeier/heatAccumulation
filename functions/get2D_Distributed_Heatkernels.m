function [heatkernel_out,q_vect_output] = get2D_Distributed_Heatkernels(heatsourceOptions,mat,timesteps)
% This function calculates the axial superposition of heat sources for various
% types of heat sources, enforcers conservation of energy, and outputs a heat kernel
% corresponding to the input specification (type, energy distribtuion, size distribution,
% evaluation position, total length of the heat source, timesteps)
% sources are mirrored about the surface boundary which is assumed to be adiabatic 
% (though losses could easily be introduced)
% the output heatkernels are normalized, multiply with the energy quantity in Joule to get
% the temperature evolution
%
% subfuntion get_internals ensures conservation of energy, interpolation etc.

% parse
x_pos = heatsourceOptions.x_pos;
y_pos = heatsourceOptions.y_pos;
z_pos = heatsourceOptions.z_pos;
len = length(heatsourceOptions.l_heat);

% check user heatsource designation vs input; adjust if necessary
if any(cell2mat(heatsourceOptions.sz_heat)>0)
    if strcmpi(heatsourceOptions.name,'point')
        error('Heatsource set to "point" but heat source sizes provided! Check your settings.')
    end
end

% outputs
q_vect_output = cell(len,1);

%% INTERNALS START HERE
%
% Now we can construct the heat kernel for the length(s) of the heat sources
% this part of the kernel is shared for all sums (prefactor / scalar)
% it is applied as last step before output
switch heatsourceOptions.name
    case {'point','ring'}
        % premultiplier for point source
        % premultiplier for ring-point source (Karkhin); note - identical to 'point'
        common_factor = 1.*(mat.rho*mat.cp*(4*pi*mat.kappa.*timesteps).^(3/2)).^-1;
    case 'gauss'
        % premultiplier for gaussian source
        common_factor = 1.*(pi*mat.rho*mat.cp*sqrt(pi*mat.kappa.*timesteps)).^-1;       
end

% calculate heatsource positions, sizes, and energy distribution,
% dependent on the number of heatsources / distance between heatsources
[heatsource_pos,sz_heat_interp,q_vect_interp] = get_internals(heatsourceOptions);

% preallocate: out_array contains the final heatkernels for the corresponding heat source lengths
heatkernel_out = zeros(len,length(timesteps));
consoleoutput = false; % by default no console output
for i = 1:len % for all provided heat source lengths
    % console output moved to here in case of single param
    if i == 1 && len > 1
        fprintf('[z,r] = %4.1f / %3.1f µm: ',z_pos*1e6,sqrt(x_pos^2+y_pos^2)*1e6)
    end
    % heatkernel calc - variant 3 - fully vectorized / broadcasted
    switch heatsourceOptions.name
        case 'point'
            % calculate mirror sources for point sources
            heatkernel = q_vect_interp{i}.'.*...
                exp(-((x_pos-0).^2+(y_pos-0).^2+(z_pos-heatsource_pos{i}.').^2)./...
                (4*mat.kappa.*timesteps));
        case 'gauss'
            % calculate mirror sources for gaussian sources
            heatkernel = q_vect_interp{i}.'.*(8*mat.kappa.*timesteps+sz_heat_interp{i}.'.^2).^-1.* ...
                exp((4*mat.kappa.*timesteps).^-1.* ...
                (((x_pos-0).^2+(y_pos-0).^2).*((sz_heat_interp{i}.'.^2.*(8*mat.kappa.*timesteps+sz_heat_interp{i}.'.^2).^-1)-1) ...
                -(z_pos-heatsource_pos{i}.').^2));
        case 'ring'
            % src: [Engineering Materials] Victor A. Karkhin - Thermal Processes in Welding (2019, Springer Singapore) pp 231-232
            % src: ring_heatsource_formula.pdf; besseli(0) -> modified bessel function first kind 0th order
            % besselfun result must only be updated if a) first loop or b) change of relevant arguments
            % importantly this is not the case for changes in the axial (z) coordinate (length of heat source)
            % this consideration can save a lot of computation time if true
            % below assumes/requires x_pos/y_pos/z_pos scalar (1,1) as is currently the case per optsheatsource
            % 
            if i == 1 || ~isequal(sz_heat_interp{i},sz_heat_interp{i-1})
                besselpart = besseli(0,sqrt((x_pos-0).^2+(y_pos-0).^2).*sz_heat_interp{i}.'./(2*mat.kappa.*timesteps));
            end
            heatkernel = besselpart.*q_vect_interp{i}.'.*...
                exp(-(sz_heat_interp{i}.'.^2+((x_pos-0).^2+(y_pos-0).^2+(z_pos-heatsource_pos{i}.').^2))./(4*mat.kappa.*timesteps));
    end
    % now sum all sources to generate the combined heatkernel
    heatkernel_sum = sum(heatkernel,1);
    % write to array
    heatkernel_out(i,:) = heatkernel_sum;
    % save current q_vect (for output we remove the internal mirror sources)
    if z_pos == 0
        q_vect_output{i} = q_vect_interp{i};
    else % remove mirror sources from output
        q_vect_output{i} = q_vect_interp{i}(1+length(q_vect_interp{i})/2:end);
    end
    % progress output if output is an array every 10% finished calc
    if len > 1 % todo: output ONLY if some time has elapsed
        if mod(i,floor(len/10)) == 0
            consoleoutput = true;
            fprintf('%1.0f%% - ',100*i/len);
        end
    end
end
% finally scale the kernel with its shared prefactor / scalar
heatkernel_out = common_factor.*heatkernel_out;
% add newline if console output was activated
if (consoleoutput == true) && (i ~= len)
    fprintf('\n'); % newline
end

end

function [heatsource_pos_out,sz_heat_out,q_vect_out] = get_internals(heatsourceOptions)

% get settings in temp vars
z_pos = heatsourceOptions.z_pos;
n_heatsources = heatsourceOptions.n_heatsources;
delta_l_heat = heatsourceOptions.delta_l_heat;
% prep output
len = length(heatsourceOptions.q_vect);
sz_heat_out = cell(len,1);
q_vect_out = cell(len,1);
heatsource_pos_out = cell(len,1);

for i = 1:len
    % get current energy distrib, length of heatsource, heatsource sizes
    q_vect = heatsourceOptions.q_vect{i};
    l_heat = heatsourceOptions.l_heat{i};
    sz_heat = heatsourceOptions.sz_heat{i};
    % get stepping (distance between point sources) from lheat
    if length(l_heat) == 1
        % when l_heat is scalar then assume linear spacing between heat sources
        source_spacing = 'calculated';
        if isnan(delta_l_heat)
            l_heat_stepping = l_heat/(n_heatsources-1); % length steps for the sum
        else
            n_heatsources = ceil(1+l_heat/delta_l_heat); % determine number of sources from heat stepping
            n_heatsources(n_heatsources < 10) = 10; % limit min amount
            if n_heatsources > 10000
                disp('Enforcing limit: No more than 10000 Heatsources')
                n_heatsources(n_heatsources > 10000) = 10000; % limit max amount
            end
            l_heat_stepping = l_heat/(n_heatsources-1); % calculate true l_heat_stepping
        end
    elseif length(l_heat) > 1
        % when l_heat is vector assume the vector describes the actual
        % positions of the heat sources
        source_spacing = 'predefined';
        n_heatsources = length(l_heat);
    end
    % generate indexes of the sum
    % special case: if evaluated at the surface then the mirror sources can
    % simply be accounted for by scalar 2 
    % for any other z_position we must evaluate the complete sum
    if z_pos == 0
        sum_idx = 0:1:(n_heatsources-1);
        scale = 2; % symmetry
    else
        sum_idx = [-(n_heatsources-1):1:0, 0:1:(n_heatsources-1)];
        scale = 1;
    end
    % calculate heat source positions
    switch source_spacing
        case 'calculated'
            heatsource_pos = l_heat_stepping.*sum_idx;
        case 'predefined'
            if z_pos == 0
                heatsource_pos = l_heat;
            elseif z_pos ~= 0
                heatsource_pos = [-fliplr(l_heat), l_heat];
            end
    end
    % parse q_vect (energy distribution along the distributed points)
    if ismember(class(q_vect),{'string','char'})
        switch lower(q_vect)
            case 'constant'
                q_vect = (1/n_heatsources)*ones(1,n_heatsources);
            otherwise
                warning('Heatkernel calc can only parse "constant" q_vect. BAD USER DEFINITION.')
                error(strcat(mfilename,'.m:',' q_vect case',32,'"',q_vect,'"',32,'unknown'));
        end
    elseif ~isvector(q_vect)
        error('q_vect must either be a known string or a vector')
    end
    % if the input q_vector does not satisfy the requirements we interpolate (ideally it is length (1,n_heatsources))
    if length(q_vect) ~= n_heatsources
        q_vect = interp1(1:length(q_vect),q_vect,linspace(1,length(q_vect),n_heatsources));
    end
    % if the input sz_heat does not satisfy the requirements we interpolate (ideally it is length (1,n_heatsources))
    if length(sz_heat) ~= n_heatsources
        if isscalar(sz_heat)
            % just for clarity: scalar is allowed and reduces complexity in heatkernel calc
        else
            sz_heat = interp1(1:length(sz_heat),sz_heat,linspace(1,length(sz_heat),n_heatsources));
        end
    end
    % normalize q_vect if it does not sum to 1 or 2 dependent on case (does not conserve energy)
    if abs(sum(q_vect)-scale) > 5*eps
        q_vect = scale.*q_vect./sum(q_vect);
    end
    % for any other position than the surface we need to evaluate the mirror
    % sources explicitly (as opposed to the special case where a scalar 2 suffices)
    % for these cases we need to construct q vect as following
    if z_pos ~= 0
        q_vect = [fliplr(q_vect), q_vect]; %#ok<AGROW>
        if ~isscalar(sz_heat)
            sz_heat = [fliplr(sz_heat),sz_heat]; %#ok<AGROW>
        end
    end
    % if a row of q_vect is all zeros, division by its sum yields NaN resulting
    % in confusing errors when applying the FFT
    if any(isnan(q_vect))
        error('Something is wrong with q_vect; NaNs have been produced')
    end
    % write output
    sz_heat_out{i} = sz_heat;
    q_vect_out{i} = q_vect;
    heatsource_pos_out{i} = heatsource_pos;
end
end