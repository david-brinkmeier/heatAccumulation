function plot3(obj)
if length(obj.sim.z_eval_depths) > 1
    % do interp
    obj.EnableInterp = 1;
    % generate figure
    fig1 = genORselectfigbyname('Temperature over depth',obj.UUID); %clf
    fig1.Position = [204 90 600 900]; fig1.Color = [1 1 1];
    % get drilldepths/hslens
    hslens = cellfun(@max,obj.sim.Heatsource2Dlen); % accounts for weber edgecase
    % select upper plot
    ax2 = subplot(2,1,1); ax2.LineWidth = 1;
    plot([1 obj.sim.TotalNoPulses],1e6*[obj.drilldepth_max obj.drilldepth_max],'--k','DisplayName','Förster')
    hold on
    plot(1:obj.sim.TotalNoPulses,1e6*hslens,'-k','DisplayName','Weber')
    ylim([0 round(1.05*1e6*obj.drilldepth_max,-1)])
    title({strcat('P_{avg,max}=',32,num2str(round(obj.MeanmaxavgLaserPower{2},2)),'W,',32,'Ep_{max}=',32,num2str(round(obj.laser.Ep*1e6)),'µJ,',32,...
        'frep_{max}=',32,num2str(max(obj.sim.reprates{1}*1e-3)),'kHz,',32,'w0 =',32,num2str(obj.laser.w0*1e6),'µm'),...
        strcat('P_{avg}=',32,num2str(round(obj.MeanmaxavgLaserPower{1},2)),'W,',32,'Q_{in,tot}=',32,num2str(round(obj.Qtotal_in,2)),'J,',32,...
        '"drill_{efficiency}"=',32,num2str(round(obj.MillimeterPerJoule,2)),32,'mm/J')})
    xlabel('Number of pulses','FontSize',12);
    ylabel('Depth in µm','FontSize',12);
    yyaxis right
    plot(1:obj.sim.TotalNoPulses,obj.gouffe,'--r','DisplayName','Gouffé'), hold on
    plot(1:obj.sim.TotalNoPulses,obj.sim.EnergyRampFull{1},'-r','DisplayName','Ep')
    plot(1:obj.sim.TotalNoPulses,obj.sim.reprates{1}./max(obj.sim.reprates{1}(:)),'-b','DisplayName','frep')
    plot(1:obj.sim.TotalNoPulses,obj.AvgLaserPower./obj.MeanmaxavgLaserPower{2},'-.','Color','#0072BD','DisplayName','Pavg')
    ylabel('Gouffé / Energy / Power / Reprate [0-1]')
    legend('Location','Southeast','NumColumns',2,'AutoUpdate','off')
    xlim([1 obj.sim.TotalNoPulses]),
    ylim([min([obj.gouffe(:); obj.sim.EnergyRampFull{1}(:)]), ...
        round(max([obj.gouffe(:); obj.sim.EnergyRampFull{1}(:)])*1.05,2)])
    % make y axes black
    ax2.YAxis(1).Color = 'k';
    ax2.YAxis(2).Color = 'k';
    % overlay
    yyaxis left
    plot1.handle = scatter(1,1e6*hslens(1),60,'s','filled');
    plot1.handle.XData = 1; plot1.handle.YData = 1e6*hslens(1);
    plot1.ydata = 1e6*hslens; % drill depth for upperplot
    % select lower plot
    plot2.axis = subplot(2,1,2); plot2.axis.LineWidth = 1;
    plot2.axis.Position = [0.13 0.15  0.77 0.35];
    plot2.tempresults.handle = plot(plot2.axis,obj.sim.z_eval_depths.*1e6,obj.TempResults_z{1,1}(:,1),'r');
    plot2.tempresults.ydata = obj.TempResults_z(1,:);
    hold(plot2.axis,'on')
    % plot max temp envelope
    plot2.tempmax.handle = plot(plot2.axis,obj.sim.z_eval_depths.*1e6,obj.MaxTemp_z(:,1),'--k');
    plot2.tempmax.ydata = obj.MaxTemp_z;
    % draw interpolated q_vects
    inith2 = ones(1,obj.interpVal); % init all ones
    plot2.qvect.handle = scatter(plot2.axis,inith2,inith2,200*inith2,'or','filled'); hold(plot2.axis,'off')
    plot2.qvect.xdata = 1e6.*obj.heat2Dlen_interp; % heat source length for q vect
    plot2.qvect.sz = 200*obj.q_vect_interp./max(obj.q_vect_interp(:)); % q vect size increase
    % set label
    ylabel('deltaT [K]'), xlabel('depth in µm');
    xlim([0 1e6*obj.sim.z_eval_depths(end)])
    % parse this data to slider listener
    titledata = {obj.sim.timesteps{1}.*1e3,obj.sim.EnergyRampFull{1}.*obj.laser.Ep.*1e6,...
        obj.sim.reprates{1}.*1e-3,obj.sim.radial_position.*1e6}; % time(n_pulse)
    
    % time slider | listener executes upon moving the slider instead of after setting the slider
    slider_time = uicontrol('Parent',fig1,'Style','slider','Position',[111,54,429,23],...
        'value',1,'min',1,'max',obj.sim.TotalNoPulses,'BackgroundColor',[1 1 1]);
    [~] = uicontrol('Parent',fig1,'Style','text','Position',[46,50,50,23],...
        'String','Process','BackgroundColor',[1 1 1]);
    % radius slider | listener executes upon moving the slider instead of after setting the slider
    slider_radial = uicontrol('Parent',fig1,'Style','slider','Position',[111,24,429,23],...
        'value',1,'min',1,'max',length(obj.sim.radial_position),'BackgroundColor',[1 1 1]);
    slider_radial_label = uicontrol('Parent',fig1,'Style','text','Position',[38,20,60,23],...
        'String','Radial pos.','BackgroundColor',[1 1 1]);
    % add listeners
    [~] = addlistener(slider_time, 'Value', 'PostSet', ...
        @(hObject,eventData) updateplot3(plot1,plot2,titledata,slider_time.Value,slider_radial.Value));
    [~] = addlistener(slider_radial, 'Value', 'PostSet', ...
        @(hObject,eventData) updateplot3(plot1,plot2,titledata,slider_time.Value,slider_radial.Value));
    if length(obj.sim.radial_position) == 1
        slider_radial.Visible = 0; slider_radial_label.Visible = 0;
    end
else
    warning('Plot3 only available if more than one z-Position is evaluated')
end
end

function updateplot3(plot1,plot2,titledata,index_t,index_r)
% plot 1 stores fighandle and data for stuff contained in upper plot
% plot 2 stores fighandles and data for stuff contained in lower plot
% titledata stores stuff for the title of the lower plot
% index_t and index_r are the slider values
%
% parse slider vals
index_t = round(index_t);
if isempty(index_r)
    index_r = 1;
else
    index_r = round(index_r);
end
% save ylims
ylim_bak = plot2.axis.YLim;
% generate and set title
titlestr = sprintf('Pulse %7.0f, Time %4.1f ms, Ep %4.1f µJ, frep %4.2f kHz, r = %3.1f µm', ...
    index_t,titledata{1}(index_t),titledata{2}(index_t),titledata{3}(index_t),titledata{4}(index_r));
title(plot2.axis,titlestr)
% for the depth in upper plot
set(plot1.handle,'XData',index_t);
set(plot1.handle,'YData',plot1.ydata(index_t));
% for the temperature over depth
set(plot2.tempresults.handle,'YData',plot2.tempresults.ydata{1,index_r}(:,index_t));
% for the max temp over depth
set(plot2.tempmax.handle,'YData',plot2.tempmax.ydata(:,index_r));
% for the q_vects
set(plot2.qvect.handle,'XData',plot2.qvect.xdata(index_t,:))
set(plot2.qvect.handle,'SizeData',plot2.qvect.sz(index_t,:));
% enforce backup ylims
plot2.axis.YLim = ylim_bak;
end