clear global
clearvars
clc
close all

%% Lade Funktionen und Pfade
addpath(genpath('functions'))
addpath(genpath('classes'))
addpath(genpath('results'))

%% Material Specification - required for all calculations

% init material
mat = materialprops('crni_1030nm');
mat.fresnel % generate fresnel absorption plot for material
% initialize laser properties
laser = initlaser('trumicro');
% initialize simulation settings
sim = simopts('init');

% Setup - Simulating a single parameter
laser.Ep = 750e-6;
laser.frep = 50e3;
laser.w0 = (60/2)*1e-6;

% simulation setup
sim.oversampling    = {1,1};                        % enable if you want to resolve Tempearature in between pulses, recommed integer >= 5
sim.TotalNoPulses   = 2^10;                         % for numerical efficiency this is ideally (2^n)/2 where n = positive integer, e.g. 1,2,3...
sim.bins            = {2^8,'static'};               % divide into [bins] chunks with constant heatkernel; static or variable allowed
sim.PulseRange      = [1 sim.TotalNoPulses];        % vector [firstpulse, lastpulse] start pulsing at firstpulse, stop at lastpulse
sim.EnergyRamp      = ones(1,sim.TotalNoPulses);    % constant pulse energy
sim.Heatsource      = 'point';                      % implemented are 'point', 'gauss' and 'ring'; heatsourcesize sets gaussian radius / ring radius respectively
sim.debug           = false;                        % adds some plots for drill depth progress and heat accumulation
sim.precision       = {'fp32','fp32'};              % 1st applies to result of storage, 2nd applies to some calculations (trade speed vs accuracy)
sim = sim.update_timesteps(sim,laser.frep,'constant'); % update timesteps for heataccu calculation according to laser settings [call after TotalNoPulses etc.]

% radial positions / grid for evaluation; 0 is ON-AXIS; everything in SI units
% sim.radial_position = 0;
% sim.radial_position = 20e-6;
% sim.radial_position = [0, 5e-6];
% sim.radial_position = 1e-6*linspace(0,100,30);
sim.radial_position = 1e-6*[0,3,6,9,12,linspace(15,100,9)];

% SETTINGS SPECIFIC TO DETERMINATION OF Z-POSITIONS TO EVALUATE (APPLIES ONLY TO "3DLine")
% sim.Depths_Eval = {0e-6,[]}; % Select a depth at which to evaluate
% sim.Depths_Eval = {20e-6,[]}; % Select a depth at which to evaluate
% sim.Depths_Eval = {[0, 20e-6],[]}; % Set explicit Z-Positions to evaluate
% sim.Depths_Eval = {linspace(0,100,10).*1e-6,[]}; % Set explicit Z-Positions to evaluate
sim.Depths_Eval = {'determine',20}; % Automatically determine [#] depths/positions to evaluate

% SETTINGS SPECIFIC TO ENERGY RAMPS
% sim.EnergyRamp = [linspace(0.01,1,sim.TotalNoPulses/2) linspace(1,1,sim.TotalNoPulses/2)];
% sim.EnergyRamp = 0.8+0.2*sin(linspace(0,6*2*pi,sim.TotalNoPulses)); % periodic
% sim.EnergyRamp = [movmean(1.5.*rand(1,sim.TotalNoPulses/2),20) linspace(1,1,sim.TotalNoPulses/2)];
% sim.EnergyRamp = sigmf(linspace(0.5,9,sim.TotalNoPulses),[1.25 3.5]);
% rampvect = [ones(1,2^3), zeros(1,2^3)]; sim.EnergyRamp = repmat(rampvect,[1, sim.TotalNoPulses/length(rampvect)]);

% SETTINGS SPECIFIC TO CHIRPED HEAT ACCUMULATION
% chirp = laser.frep*ones(1,sim.TotalNoPulses); % override constant reprate force explicit calc
% chirp = [linspace(300e3/9,150e3,sim.TotalNoPulses/2) 150e3*ones(1,sim.TotalNoPulses/2)];
% chirp = [linspace(300e3/9,150e3,sim.TotalNoPulses/1.5) 150e3*ones(1,ceil(sim.TotalNoPulses/3))];
% chirp = [linspace(10e3,25e3,sim.TotalNoPulses/2) 25e3*ones(1,sim.TotalNoPulses/2)];
% chirp = 1e3.*[round((150-35).*rand(1,sim.TotalNoPulses/2),0)+35, 150*ones(1,sim.TotalNoPulses/2)]; 
% sim = sim.update_timesteps(sim,chirp,'variable');
% sim = sim.update_timesteps(sim,chirp,'variable-intdivider',laser.frep); % only allow integer dividers applied to base/reference pulse repetition rate

% SETTINGS SPECIFIC TO RELATIVE ENERGY DISTRIBUTION ALONG LINE (APPLIES ONLY TO "3DLine")
sim.Q_distrib_2D = {'drilldepthmodel'}; % corresponding to linear fluence increase hth -> htip
% sim.Q_distrib_2D = {'constant'};
% sim.Q_distrib_2D = {[zeros(1,9),1]};
% sim.Q_distrib_2D = {linspace(1,4,40)};
% sim.Q_distrib_2D = {linspace(1,4,20).*repmat([0,0,0,0,1],[1,4])};
% sim.Q_distrib_2D = {linspace(eps,2,30)};
% sim.Q_distrib_2D = {[2,zeros(1,38),1]}; % surface and focus heating
% sim.Q_distrib_2D = {[0.5,zeros(1,11),1,zeros(1,13),2,zeros(1,12),4]}; % surface and focus heating
% sim.Q_distrib_2D = {rand(sim.TotalNoPulses,30)};
% sim.Q_distrib_2D = {linspace(0.5,1,sim.TotalNoPulses).'.*linspace(0,1,50)};
% sim.Q_distrib_2D = repmat({linspace(1,2,30)},[sim.TotalNoPulses 1]);

% explicity Heat source length for/after each Pulses (this disables drilling progress calculation
% sim.Heatsource2Dlen = {300e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
% sim.Heatsource2Dlen = num2cell(1e-6.*(100+70.*sin(linspace(0,8*2*pi,sim.TotalNoPulses))).',2); % periodic changes

% SETTINGS SPECIFIC TO GENERATION OF LINE SOURCE FROM 3D POINT SOURCES 
% sim.HeatsourceSettings.n_heatsources = 25;
sim.HeatsourceSettings.delta_l_heat = 15e-6;
% sim.HeatsourceSettings.delta_l_heat = 0.2*get_lthermal(mat,min(sim.timesteps{1})).lthermal; % make it depend on thermal difflen

% set beam radius which affecs heat source geometry (default == 0 / Point Source)
% sim.Heatsource2Dsize = {0e-6}; % same as default value
% sim.Heatsource2Dsize = {1e-6.*[0.1,0.1]}; % compare with point source
% sim.Heatsource2Dsize = {1e-6.*[(40/2),(40/2)]}; % const w0
% sim.Heatsource2Dsize = {1e-6.*10}; % const w0
% sim.Heatsource2Dsize = num2cell(1e-6.*(50+30.*sin(linspace(0,5*2*pi,sim.TotalNoPulses))).',2); % periodic changes
% sim.Heatsource2Dsize = {1e-6.*linspace(42.9/2,sqrt(2)*(42.9/2),40)}; % w(0) focus on workpiece, after w(zR) = sqrt(2)*w(0)
% sim.Heatsource2Dsize = {1e-6.*linspace(100,30,20)}; % 200µm on surface, 10 µm in focus
% sim.Heatsource2Dsize = {1e-6.*linspace(10,100,40).*randi([0 1],1,40)}; % 200µm on surface, 10 µm in focus
% sim.Heatsource2Dsize = num2cell(1e-6*linspace(1,50,10).*linspace(1,2,sim.TotalNoPulses).',2); % defines Gaussian / Ringsource "caustic" for / @ each pulse
% sim.Heatsource2Dsize = {1e-6*linspace(1,50,10).*linspace(1,2,sim.TotalNoPulses).'}; % defines Gaussian / Ringsource "caustic" for / @ each pulse

% calculate
% profile on
ttot = tic;
results = heatacc(laser,mat,sim);
fprintf('Total calc time: %5.1fs, heat conduction: %5.1fs\n',toc(ttot),results.totalsimtime);
% profile viewer

%% variants of plot export
% close all

results.plotall;
% results.plot1;
% results.plot2;
% results.plot3;
% results.plot4;
% plotcfg = results.plot4; % start plot 4 and save export config
% plotcfg = results.plotall; % start all plots and save config of plot 4
% results.plot4(plotcfg); % pass export config to other results as initial (reproduce exact video output)


%% TESTCASES

%{
mat = materialprops('crni_1030nm');
% initialize laser properties
laser = initlaser('trumicro');
laser.Ep = 1500e-6;
laser.frep = 25e3;
laser.w0 = (107/2)*1e-6;
% initialize simulation settings
sim = simopts('init');

% simulation setup
sim.oversampling    = {1,1};                        % enable if you want to resolve Tempearature in between pulses, recommed integer >= 5
sim.TotalNoPulses   = 2^11;                         % for numerical efficiency this is ideally (2^n)/2 where n = positive integer, e.g. 1,2,3...
sim.EnergyRamp      = ones(1,sim.TotalNoPulses);    % constant pulse energy
sim.bins            = {2^7,'variable'};             % divide into [bins] chunks with constant heatkernel; static or variable allowed
sim.PulseRange      = [1 sim.TotalNoPulses-100];     % vector [firstpulse, lastpulse] start pulsing at firstpulse, stop at lastpulse
sim.Heatsource      = 'point';                     % implemented are 'point', 'gauss' and 'ring'; heatsourcesize sets gaussian radius / ring radius respectively
sim.debug           = false;                        % adds some plots for drill depth progress and heat accumulation
sim.precision       = {'fp32','fp32'};              % 1st applies to result of storage, 2nd applies to some calculations (trade speed vs accuracy)
sim = sim.update_timesteps(sim,laser.frep,'constant'); % update timesteps for heataccu calculation according to laser settings [call after TotalNoPulses etc.]

% z,r positions
sim.radial_position = 1e-6*linspace(0,110,5);
sim.Depths_Eval = {'determine',20}; % Automatically determine depths 
% number of heatsources
sim.HeatsourceSettings.n_heatsources = 20;

% init
tcomplete = tic;
results = cell({}); 
idx = 0;

% starting here
idx = idx+1;
sim.Heatsource = 'point';
sim.Q_distrib_2D = {'constant'};
sim.Heatsource2Dsize = {0e-6};
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
results{idx,1} = heatacc(laser,mat,sim);
%
idx = idx+1;
sim.Heatsource = 'gauss';
sim.Q_distrib_2D = {'constant'};
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
sim.Heatsource2Dsize = {1e-6.*0.1}; % compare with point source
results{idx,1} = heatacc(laser,mat,sim);
%
idx = idx+1;
sim.Heatsource = 'ring';
sim.Q_distrib_2D = {'constant'};
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
sim.Heatsource2Dsize = {1e-6.*0.1}; % compare with point source
results{idx,1} = heatacc(laser,mat,sim);
%
idx = idx+1;
sim.Heatsource = 'point';
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
sim.Q_distrib_2D = num2cell(linspace(0.5,1,sim.TotalNoPulses).'.*linspace(0,1,50),2);
sim.Heatsource2Dsize = {0}; % defines Gaussian / Ringsource "caustic" for / @ each pulse
results{idx,1} = heatacc(laser,mat,sim);
%
idx = idx+1;
sim.Heatsource = 'gauss';
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
sim.Q_distrib_2D = num2cell(linspace(0.5,1,sim.TotalNoPulses).'.*linspace(0,1,50),2);
sim.Heatsource2Dsize = num2cell(1e-6*linspace(1,50,10).*linspace(1,2,sim.TotalNoPulses).',2); % defines Gaussian / Ringsource "caustic" for / @ each pulse
results{idx,1} = heatacc(laser,mat,sim);
%
idx = idx+1;
sim.Heatsource = 'ring';
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
sim.Q_distrib_2D = num2cell(linspace(0.5,1,sim.TotalNoPulses).'.*linspace(0,1,50),2);
sim.Heatsource2Dsize = num2cell(1e-6*linspace(1,50,10).*linspace(1,2,sim.TotalNoPulses).',2); % defines Gaussian / Ringsource "caustic" for / @ each pulse
results{idx,1} = heatacc(laser,mat,sim);
%
% % now switching to variable length heat sources
%
sim.HeatsourceSettings.delta_l_heat = 15e-6;
% second run here
idx = idx+1;
sim.Heatsource = 'point';
sim.Q_distrib_2D = {'constant'};
sim.Heatsource2Dsize = {0e-6};
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
results{idx,1} = heatacc(laser,mat,sim);
%
idx = idx+1;
sim.Heatsource = 'gauss';
sim.Q_distrib_2D = {'constant'};
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
sim.Heatsource2Dsize = {1e-6.*0.1}; % compare with point source
results{idx,1} = heatacc(laser,mat,sim);
%
idx = idx+1;
sim.Heatsource = 'ring';
sim.Q_distrib_2D = {'constant'};
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
sim.Heatsource2Dsize = {1e-6.*0.1}; % compare with point source
results{idx,1} = heatacc(laser,mat,sim);
%
idx = idx+1;
sim.Heatsource = 'point';
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
sim.Q_distrib_2D = num2cell(linspace(0.5,1,sim.TotalNoPulses).'.*linspace(0,1,50),2);
sim.Heatsource2Dsize = {0}; % defines Gaussian / Ringsource "caustic" for / @ each pulse
results{idx,1} = heatacc(laser,mat,sim);
%
idx = idx+1;
sim.Heatsource = 'gauss';
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
sim.Q_distrib_2D = num2cell(linspace(0.5,1,sim.TotalNoPulses).'.*linspace(0,1,50),2);
sim.Heatsource2Dsize = num2cell(1e-6*linspace(1,50,10).*linspace(1,2,sim.TotalNoPulses).',2); % defines Gaussian / Ringsource "caustic" for / @ each pulse
results{idx,1} = heatacc(laser,mat,sim);
%
idx = idx+1;
sim.Heatsource = 'ring';
sim.Heatsource2Dlen = {200e-6.*ones(sim.TotalNoPulses,1)}; % fixed heat source length
sim.Q_distrib_2D = num2cell(linspace(0.5,1,sim.TotalNoPulses).'.*linspace(0,1,50),2);
sim.Heatsource2Dsize = num2cell(1e-6*linspace(1,50,10).*linspace(1,2,sim.TotalNoPulses).',2); % defines Gaussian / Ringsource "caustic" for / @ each pulse
results{idx,1} = heatacc(laser,mat,sim);
% write complete calculation time
results{idx+1,1} = toc(tcomplete);
%}
