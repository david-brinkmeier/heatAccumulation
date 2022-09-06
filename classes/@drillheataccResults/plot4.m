function plotcfg = plot4(obj,plotcfg)
%#ok<*TRYNC>

skip = false;
if any(size(obj.TempResults_zrt,1:2) == 1)
    warning('Plot4 requires >1 r,z positions, respectively. Skipping.')
    skip = true;
end

if exist('plotcfg','var') && isstruct(plotcfg)
    plotdefprovided = true;
else
    plotdefprovided = false;
end

plottype = [];
if ~skip
plottype = lower(questdlg('Specify plot type for plot 4.', ...
    'Plot type selector', ...
    'Contour','Surface','Surface'));
end

if ~isempty(plottype) && ~skip % abort selector
    % prepare fig and axes
    data.fig = genORselectfigbyname('T(r,z) - interactive / setup',obj.UUID); clf(data.fig); set(data.fig,'units','pixel')
    data.fig.Color = [1 1 1]; data.fig.Position = [810 600 540 720];
    data.ax = axes(data.fig,'Position',[0.1,0.05,0.72,0.815],...
        'OuterPosition',[0.03,0.03,0.95,0.95],'LineWidth',1);
    data.fontsz = [15,12]; % fonsizes
    % prepare data
    data.UUID = obj.UUID;
    data.img.zdata = obj.TempResults_zrt;
    startimg = [flip(data.img.zdata(:,2:end,1),2), data.img.zdata(:,:,1)];
    data.plottypes = {'matlab.graphics.chart.primitive.Contour',...
        'matlab.graphics.chart.primitive.Surface'};
    data.titlestr = {'t = %1.1f ms, E_{res} = %3.1f µJ, E_{in} = %3.1f µJ\nf_{rep} = %3.1f kHz, P_{in} = %3.1f W, l_{hs} = %3.1f µm',...
        't = %1.1f ms, E_{res} = %3.1f µJ\n E_{in} = %3.1f µJ, f_{rep} = %3.1f kHz\nP_{in} = %3.1f W, l_{hs} = %3.1f µm',1};
    data.titledata = {obj.sim.timesteps{1}.*1e3,obj.E_res.*1e6,...
        obj.E_in.*1e6,obj.sim.reprates{1}.*1e-3,...
        obj.E_in.*obj.sim.reprates{1}};
    data.img.hslen = cellfun(@(x) x(end)*1e6,obj.sim.Heatsource2Dlen); % heatsource lengths
    data.effective_pulse = cumsum(logical(obj.sim.EnergyRampFull{1}.*obj.sim.HeatInputs{1})); % effective number of applied pulses
    data.img.zdatamax = obj.MaxTemp_z;
    data.img.hslenmax = data.img.hslen(end);
    data.img.radiusvect = round(1e6.*[-flip(obj.sim.radial_position(2:end)), obj.sim.radial_position],1); % radial positions
    data.img.depthvect = round(obj.sim.z_eval_depths.*1e6,1); % depths
    data.stretch_horizontal = 1; % daspect stretch
    data.numofxticks = 5; % number of xticks
    data.xlims = [min(data.img.radiusvect),max(data.img.radiusvect)];
    data.ylims = [min(data.img.depthvect),max(data.img.depthvect)];
    data.tcrit = obj.mat.Tcritical;
    data.templims = [0, data.tcrit];
    
    % make not whether it's drilling process
    data.isdrilling = obj.sim.isdrilling;
    
    % prepare annotation
    data.annotation = annotation('textbox', [0, 0, 0.5, 0],...
        'Fontsize',data.fontsz(1),'Fontweight','normal','Color','#cb38cb',...
        'EdgeColor','none','string',sprintf('%4.0f',data.effective_pulse(1)));
    
    % draw first frame / prepare colormap if contourf
    switch plottype
        case 'contour'
            data.contour.levels = 10;
            data.contour.isolines = round(linspace(data.templims(1),data.templims(2),data.contour.levels),-1);
            customcmap = interp2(jet(255),1:3,linspace(1,255,data.contour.levels-1).');
            warning('off','MATLAB:contour:ConstantData')
            [~,data.plot] = contourf(data.ax,data.img.radiusvect,data.img.depthvect,startimg,data.contour.isolines,'ShowText','off');
            colormap(data.ax,customcmap)
            warning('on','MATLAB:contour:ConstantData')
        case 'surface'
            data.plot = surf(data.ax,data.img.radiusvect,data.img.depthvect,startimg); view(2), shading interp
            data.surface.levels = 8;
            colormap(data.ax,jet(255))
    end
    
    hold on, % add one pretty line on top
    plot3(data.ax,[min(data.img.radiusvect),max(data.img.radiusvect)],[0,0],[1e7,1e7],'-k')
    hold off
    
    % setup figure / labels / fontsizes
    set(data.ax,'Ydir','reverse')
    xlim(data.ax,data.xlims)
    ylim(data.ax,data.ylims)
    
    % fix ticks
    xticksvect = round(linspace(0,data.img.radiusvect(end),data.numofxticks),0);
    xticksvect = [-flip(xticksvect(2:end)), xticksvect];
    xticks(data.ax,xticksvect);
    xticklabels(data.ax,num2cell(abs(xticksvect.'))); % no negative radii
    
    data.cbar = colorbar(data.ax);
    caxis(data.ax,[0 round(data.templims(2),-1)]); %1500+557 where 557 is melt enthalpy
    switch plottype
        case 'contour'
            set(data.cbar,'YAxisLocation','right','YTick',data.contour.isolines);
        case 'surface'
            set(data.cbar,'YAxisLocation','right');
    end
    data.labels.x = xlabel(data.ax,'radius in µm','FontSize',data.fontsz(1)+2);
    data.labels.y = ylabel(data.ax,'depth in µm','FontSize',data.fontsz(1)+2);
    data.labels.cbar = ylabel(data.cbar, 'temperature increase in K','FontSize',data.fontsz(1));
    data.cbar.FontSize = data.fontsz(2); data.labels.cbar.FontSize = data.fontsz(1);
    data.ax.XAxis.FontSize = data.fontsz(1); data.labels.x.FontSize = data.fontsz(1)+2;
    data.ax.YAxis.FontSize = data.fontsz(1); data.labels.y.FontSize = data.fontsz(1)+2;
    daspect(data.ax,[1 1 1]) % set aspect ratio, allow for stretching (2)
    
    % init patch if drilldepth
    if data.isdrilling(1)
        coord.x = 1e6.*[-data.isdrilling(2) data.isdrilling(2) 0]; coord.y = zeros(1,3);
        coord.z = 10.*repelem(max(obj.MaxTemp_z(:)),3);
        data.patch = patch(data.ax,coord.x,coord.y,coord.z,...
            'white','LineStyle','-','LineWidth',0.5,'EdgeColor','k');
    end
    
    % init frames
    data.frames.currentframe = 1;
    data.frames.actualframe = 1;
    data.frames.start = 1;
    data.frames.end = obj.sim.TotalNoPulses;
    data.frames.lastframe = obj.sim.TotalNoPulses;
    data.frames.step = 1;
    data.frames.showmax = 0; % forces on/off max temp
    
    % video export spec
    data.export.enable = false;
    data.export.fps = 30;
    data.export.codec = 'mp4';
    
    % now just overwrite data if plotcfg is provided...
    if plotdefprovided
        try data.fontsz = plotcfg.fontsz; end
        try data.stretch_horizontal = plotcfg.stretch_horizontal; end
        try data.numofxticks = plotcfg.numofxticks; end
        try data.ylims = plotcfg.ylims; end
        try data.contour = plotcfg.contour; end
        try data.frames.start = plotcfg.frames.start; end
        try data.frames.end = plotcfg.frames.end; end
        try data.frames.step = plotcfg.frames.step; end
        try data.export = plotcfg.export; end
        try data.fig.Position = plotcfg.figsz; end
        try data.templims = plotcfg.templims; end
        try data.surface.levels = plotcfg.surface.levels; end
        % dont judge me gui stuff is taking too much time
    end
    
    % IMPORTANT: store data in the figure! (retrieve: data = guidata(data.fig);
    guidata(data.fig,data)
    
    % now start gui
    try
        plotcfg = drillheataccResults.plot4gui(data.fig,obj.UUID,obj.sim.TotalNoPulses,plotdefprovided);
    end
else
    plotcfg = []; 
end

end