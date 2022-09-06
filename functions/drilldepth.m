function [depth_array_output,Abs_Fact,Gouffe_Abs,drilldepth_max,energy_distribution,ablation_crater_radius] = drilldepth(laser,mat,sim)
% Requires class initlaser and materialprops
% Rampvect is a (1,n_pulses) vector with values ranging from 0 to 1 used to
% describe energy ramps
%
% Ep: Pulse energy in Joule
% w0: Beam radius in microns
% Np: Number of pulses to calculate
% rampvect: Rampvect is a vector of length(Np) of value 0<= val<= 1 and
% may be provided to consider energy ramps
% (current implementation does not consider pulse energy dependent / varying rthh)
% debug: Debug flag activates plot
%
% Output depth_array_output is depth in SI units [m]

% Description:
% Script that calculates the drilling depth based on an analytical model
% (not expüerimentally validated for varying pulse energy etc.)
% (not valid if significant melting occurs)
% (residual energy distribution over borehole depth according to assumptions of 
% drill progress depth model is at odds with experimental observations)
% YMMV

%% Parse inputs / change to the units used in this function

Ep = laser.Ep*1e3;      % convert pulse energy to mJ
w0 = laser.w0*1e6;      % convert laser radius to µm
hV = mat.hV;            % evaporization energy in J/mm^3
lalpha = mat.pdepth*1e9; % energy absorption depth [nm]
A = mat.A0;        % initial absorption coefficient for drill model

n_pulses = sim.TotalNoPulses;                                           % Total number of pulses
rampvect = sim.EnergyRampFull{1};                                       % Vector (1,n_pulses) with values 0-1 where 1 = laser.Ep
reprates = sim.reprates{1}./max(sim.reprates{1}(:));                    % only used for plotting chirps  
avg_pwr = (sim.EnergyRampFull{1}.*sim.HeatInputs{1}.*laser.Ep)./ ...
           sim.timestepping{1};                                         % average laser power
mean_pwr = sum(sim.EnergyRampFull{1}.*sim.HeatInputs{1}.*laser.Ep)/...  % mean laser power until process end
           sim.timesteps{1}(sim.PulseRange(2));
max_pwr = max(avg_pwr);                                                 % maximum laser power
pwr_normalized = avg_pwr./max_pwr;                                      % normalized laser power
total_q_in = sum(sim.EnergyRampFull{1}.*...                             % total energy input in Joule
                 sim.HeatInputs{1}.*laser.Ep);
debug = sim.debug;                                                      % Debug flag (produces plot)

%% Error checks

if w0 == 0
    error('drill depth and 2d heat accumulation cannot be calculated with w0 = 0µm (point source)');
end

% if no energy ramp is provided simply fill a vector of correct length with ones
if isempty(rampvect)
    rampvect = ones(n_pulses,1);
end

% ablation crater rthh / geometry wildly varies with pulse energy
% this model can currently only properly handle constant pulse energies AND cooling
% cycles (i.e. no pulse energy); output warning in case of energy ramps
if max(abs(diff(rampvect))) ~= 0
    warning('This drill depth model cannot gracefully handle varying pulse energies other than ON/OFF cycles. YMMV.');
end

%% Calculate current and threshold fluence
phi_th = lalpha*10^-6*hV*10^2; % J/cm²
phi_0 = 2*Ep/(pi*w0^2)*10^5; % J/cm²

%% model by Förster2018 Estimation of the depth limit
F_z_drill = sqrt((phi_0^2-phi_th^2*log(phi_0/phi_th)^2)/(2*phi_th^2*log(phi_0/phi_th)))*w0;

warning('Change foerster depth definition to new one')
% change foerster depth to new definition                   
% H0 = 2*Havg*abs;
% eta_gouffe = 1;
% zlim_quality =    d0*sqrt(...
%                 ( ((eta_gouffe/abs).*H0).^2.*(1-Hth_A./H0).^2 - Hth_A^2*log(H0./Hth_A).^2 )...
%                 ./...
%                 (8*Hth_A^2*log(H0/Hth_A)) );

%% calculation
depth = 0.001*10^-3; % nm -> initial borehole depth;
depth_array = zeros(n_pulses,1); % depth progress per pulse
phi_tip = zeros(n_pulses,1); % tip fluences, required for energy distrib
Gouffe_Abs = zeros(n_pulses,1); % Gouffé Absorption array
Abs_Fact = zeros(n_pulses,1); % Factor Gouffé / Absorption
ablation_crater_radius = zeros(n_pulses,1); % store ablation crater

disablemsg = false;
for n = 1:n_pulses
    if rampvect(n) > 0
        rthh = (w0*10^-3*sqrt(0.5*log(2*rampvect(n)*Ep*10^-3/(pi*(w0*10^-3)^2*lalpha*10^-6*hV))))*10^3; % ablation crater radius in µm
        ablation_crater_radius(n) = rthh;
    end
    % Gouffé
    omegaG = 4*pi*(sin(0.5*atan(2*rthh/(2*depth))))^2;
    OoK = 1/(1+sqrt(1+(2*depth/(2*rthh))^2));
    A_GoK = A*(1+(1-A)*(OoK-omegaG/(2*pi)))/(A*(1-OoK)+OoK); % Gouffé Absorption
    
    % progress according to Holder et al. https://doi.org/10.1007/s00339-021-04455-3
    phi_tip(n,1) = (3*A_GoK*(rampvect(n)*Ep*1e-3))/(1e-8*pi*rthh*sqrt(rthh^2+depth^2))-2*phi_th; % fluence at tip in J/cm²
    zabl = lalpha*10^-3*log(phi_tip(n,1)/phi_th); % ablation depth at tip
    depth = depth + zabl; % current absolute drill depth
    
    % handling of special cases
    if ((rampvect(n) == 0 || any(imag([rthh omegaG A_GoK]))) || zabl < 0)&& n > 1
        depth_array(n,1) = depth_array(n-1,1);
        depth = depth_array(n,1);
        phi_tip(n,1) = 0;
        Gouffe_Abs(n,1) = Gouffe_Abs(n-1,1);
        Abs_Fact(n,1) = Abs_Fact(n-1,1);
    elseif n == 1
        depth_array(n,1) = 0.001*10^-3; % set to start drill depth
        Gouffe_Abs(n,1) = A;
        Abs_Fact(n,1) = 1;
    else
        depth_array(n,1) = depth; 
        Gouffe_Abs(n,1) = A_GoK;
        Abs_Fact(n,1) = A_GoK/A;
    end

    if ~(isreal(zabl) && zabl > 0) && (disablemsg == false) && (rampvect(n) ~= 0)
        msg = '(Initial) pulse energy too low, no drilling. Imaginary or negative drill depth result. Consider adjusting pulse energy or energy ramp.';
        % errordlg(msg);
        warning(msg);
        disablemsg = true;
    end
    
end

% currently only output max ablation crater
ablation_crater_radius = max(ablation_crater_radius)*1e-6;

% output maximum analytical depth, conert to SI-units [m]
drilldepth_max = F_z_drill*1e-6;

% convert output depth to SI-units [m]
depth_array_output = depth_array*1e-6;

% calculate depth per energy in µm/Joule
MikrometerPerJoule = (max(depth_array_output)*1e6)/total_q_in;

%% hacky implementation to calculate energy distribution along sources
% check fluence_vs_energy_drilldepth_model.m for reference

if strcmp(sim.Q_distrib_2D{:},'drilldepthmodel')
    fluence_z = @(htip,hth,depth,z) (((htip-hth)/depth).*z+hth);
    r_avg = @(r_abl,depth,z) (-r_abl/depth).*z+r_abl;
    
    energy_distribution = zeros(n,50); % init
    for i = 1:n_pulses
        % calculate only if flunce / pulse energy > 0
        if phi_tip(i,1) ~= 0
            z_vect = linspace(0,depth_array(i,1),50); % evaluated depths µm
            r_avg_distrib = r_avg(rthh,depth_array(i,1),z_vect); % average corresponding radius µm
            fluence_distrib = fluence_z(phi_tip(i,1),phi_th,depth_array(i,1),z_vect); % local fluence in J/cm²
            dz = z_vect(2)-z_vect(1); % gradient z
            dr = r_avg_distrib(2)-r_avg_distrib(1); % gradient radial
            
            % apply local fluence to illuminated cone shell / cylinder
            % delta_energy = fluence_distrib.*r_avg_distrib*2*pi*dz; % approximation: illuminated cylinder
            delta_energy = fluence_distrib.*(2*r_avg_distrib+abs(dr)).*pi.*sqrt(dr^2+(dz^2)); % illuminated cone shell area = (R+r)*pi*m, where m = sqrt((R-r)^2+h^2);
            
            % conversavtion of energy
            energy_distribution(i,:) = delta_energy./(sum(delta_energy(:))); % normalize
        elseif i == 1
            % else just write ones...it's deactivated anyway bc of pulsematrix
            energy_distribution(i,:) = ones(1,50)./50; % set to ones
        else
            % beware: when binning is used a mean of energy distributions
            % is used for the superposition. If the energy distribution wildly varies
            % from bin to bin this can produce artifacts. Energy distribution is crucial.
            energy_distribution(i,:) = energy_distribution(i-1,:); % just use last one
        end
    end
else
    energy_distribution = NaN;
end

% energy_distribution = energy_distribution+(0.6*max(energy_distribution(:)).*rand(n,50));
% warning('random energy distrib activated')
% val = ceil(n_pulses/100); clf; plot(energy_distribution(val,:)); % debug / plotting

%% plots and data
if debug == true
    genORselectfigbyname('Drill depth progress',sim.UUID);
    clf
    plot([1 n_pulses],[F_z_drill F_z_drill],'--k','DisplayName','Förster')
    hold on
    plot(1:1:n_pulses,depth_array,'-k','DisplayName','Weber')
    ylim([0 round(1.05*F_z_drill,-1)])
    title({strcat('P_{avg,max}=',32,num2str(round(max_pwr,2)),32,'W,',32,'Ep_{max}=',32,num2str(round(Ep*1e3)),32,'µJ,',32,...
                 'frep_{max}=',32,num2str(max(sim.reprates{1}*1e-3)),'kHz,',32,'w0 =',32,num2str(w0),32,'µm'),...
           strcat('P_{avg}=',32,num2str(round(mean_pwr,2)),32,'W,',32,'Q_{in,tot}=',32,num2str(round(total_q_in,2)),32,'J,',32,...
           '"drill_{efficiency}"=',32,num2str(round(MikrometerPerJoule,2)),32,'µm/J')})
    xlabel('Number of pulses','FontSize',12);
    ylabel('Depth in µm','FontSize',12);
    yyaxis right
    plot(1:1:n_pulses,Gouffe_Abs,'--r','DisplayName','Gouffé')
    hold on
    plot(1:1:n_pulses,rampvect,'-r','DisplayName','Ep')
    plot(1:1:n_pulses,reprates,'-b','DisplayName','frep')
    plot(1:1:n_pulses,pwr_normalized,'-.','Color','#0072BD','DisplayName','Pavg')
    ylabel('Gouffé / Energy / Power / Reprate [0-1]')
    legend('Location','Southeast','NumColumns',2)
    xlim([1 n_pulses]), ylim([min([Gouffe_Abs(:); rampvect(:)]), round(max([Gouffe_Abs(:); rampvect(:)])*1.05,2)])
    
    % make y axes black
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    drawnow
end

end