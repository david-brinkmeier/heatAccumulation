function plot2(obj)
% we pass pulse and time range
fig = genORselectfigbyname('Heat accumulation',obj.UUID);
fig.Color = [1 1 1]; fig.Position = [810 90 560 420];
% prepare data
plot1.ydata = obj.TempResults_z(1,:); % temp over z cell array
plot1.xrange = {1:obj.sim.TotalNoPulses,1e3.*obj.sim.timesteps{1}}; % pulses / time
plot1.labels = {'Number of pulses','Time in ms'};
plot1.zpos = obj.sim.z_eval_depths;
plot1.rpos = obj.sim.radial_position;
plot1.axis = axes('Position',[0.125 0.25 0.775 0.675]);
plot1.handle = plot(plot1.xrange{1},plot1.ydata{1}(1,:),'-k','DisplayName',obj.sim.Heatsource); hold on
% figure details / static stuff
plot([0 plot1.xrange{1}(end)],[obj.mat.Tcritical obj.mat.Tcritical],'--r','DisplayName','Tcrit'), hold off
xlim([plot1.xrange{1}(1) plot1.xrange{1}(end)])
xlabel(plot1.labels{1}), ylabel('Temperature in °C')
ax = gca; ax.LineWidth = 1;
legend('AutoUpdate','off')
title(sprintf('[z,r] = %4.1f / %3.1f µm: ',obj.sim.z_eval_depths(1)*1e6,obj.sim.radial_position(1)*1e6))
drawnow
% generate toggle button
button1 = uicontrol('style','togglebutton',...
    'units','pix',...
    'position',[10 8 60 23],...
    'fontsize',12,...
    'BackgroundColor',[1 1 1],...
    'string','switch');
% regular callback
button1.Callback = @(hObject,eventData) updateplot2_1(plot1,button1.Value);
% generate z-slider
slider_zpos = uicontrol('Parent',fig,'Style','slider','Position',[141,8,366,23],...
    'value',1,'min',1,'max',length(obj.sim.z_eval_depths),'BackgroundColor',[1 1 1]);
slider_zpos_label = uicontrol('Parent',fig,'Style','text','Position',[81,8,60,19],...
    'String','z-Position:','BackgroundColor',[1 1 1]);
slider_rpos = uicontrol('Parent',fig,'Style','slider','Position',[141,36,366,23],...
    'value',1,'min',1,'max',length(obj.sim.radial_position),'BackgroundColor',[1 1 1]);
slider_rpos_label = uicontrol('Parent',fig,'Style','text','Position',[81,36,60,19],...
    'String','r-Position:','BackgroundColor',[1 1 1]);
% listener executes upon moving the slider instead of after setting the slider
[~] = addlistener(slider_zpos, 'Value', 'PostSet',...
    @(hObject,eventData) updateplot2_2(plot1,slider_zpos.Value,slider_rpos.Value));
[~] = addlistener(slider_rpos, 'Value', 'PostSet',...
    @(hObject,eventData) updateplot2_2(plot1,slider_zpos.Value,slider_rpos.Value));
% / listener
if length(obj.sim.z_eval_depths) < 2 % then hide z slider
    slider_zpos.Visible = 0; slider_zpos_label.Visible = 0;
    slider_rpos.Position = [141,8,366,23]; % switch rz pos
    slider_rpos_label.Position = [81,8,60,19]; % switch rz pos
end
if length(obj.sim.radial_position) < 2 % then hide r slider
    slider_rpos.Visible = 0; slider_rpos_label.Visible = 0;
end
end

function updateplot2_1(plot1,buttonval)
% plot1 struct holds all handles and required data
buttonval = buttonval+1; % toggle switch is 0-1 we need 1-2
plot1.handle.XData = plot1.xrange{buttonval};
plot1.axis.XLim = ([plot1.xrange{buttonval}(1) plot1.xrange{buttonval}(end)]);
plot1.axis.XLabel.String = plot1.labels{buttonval};
end

function updateplot2_2(plot1,index_z,index_r)
% plot1 struct holds all handles and required data
% index_z and index_r are the slider values
index_z = round(index_z);
index_r = round(index_r);
% generate and set title
title(plot1.axis,sprintf('[z,r] = %4.1f / %3.1f µm: ',plot1.zpos(index_z)*1e6,plot1.rpos(index_r)*1e6))
% for the depth in upper plot
set(plot1.handle,'YData',plot1.ydata{index_r}(index_z,:));
end
