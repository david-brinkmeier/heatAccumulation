function output = plot4gui(referencefig,uuid,lastframe,plotdefprovided)
%#ok<*NASGU>

% init gui
figname = strcat('T(r,z) controls -',32,uuid);
guifig = figure('name',figname,'color','w');
set(guifig,'MenuBar','none'); set(guifig,'ToolBar','none');
guifig.Position = [160,160,800,300]; guifig.Resize = 0;

% read data from fig (required if plotdefprovided requires overwriting initial values
data = guidata(referencefig);

% store some local gui variables in the gui fig
% https://de.mathworks.com/help/matlab/ref/setappdata.html
setappdata(guifig,'playing',0);
setappdata(guifig,'forcestop',0);
setappdata(guifig,'framestep',1);
setappdata(guifig,'lastframe',lastframe);
setappdata(guifig,'framerange_export',[1,lastframe]);
setappdata(guifig,'ylimsinit',[min(data.img.depthvect),max(data.img.depthvect)]);
setappdata(guifig,'currentylims',[min(data.img.depthvect),max(data.img.depthvect)]);
setappdata(guifig,'stretchval',1); % data aspect ratio horizontal
setappdata(guifig,'vidfps','30');
setappdata(guifig,'vidcodec','mp4');
setappdata(guifig,'mat_tcrit',data.tcrit);
setappdata(guifig,'currentclims',[0,data.tcrit]);

if plotdefprovided
    setappdata(guifig,'framerange_export',[data.frames.start,data.frames.end]);
    setappdata(guifig,'framestep',data.frames.step);
    setappdata(guifig,'currentylims',[data.ylims(1),data.ylims(2)]);
    setappdata(guifig,'stretchval',data.stretch_horizontal);
    setappdata(guifig,'vidfps',num2str(data.export.fps));
    setappdata(guifig,'vidcodec',data.export.codec);
    setappdata(guifig,'currentclims',data.templims);
end

% main panel
mainpanel = uipanel('Title','main.settings','FontSize',12,...
             'BackgroundColor','white','Position',[.025 .025 .95 .4]);
slider_time = uicontrol('Parent',mainpanel,'Style','slider',...
    'Units','normalized','Position',[.025,.5,0.95,.4],...
    'Tag','slider_time',...
    'value',1,'min',1,'max',lastframe,...
    'BackgroundColor',[1 1 1]);
playbutton = uicontrol('Parent',mainpanel,'Style','togglebutton',...
    'Units','normalized','Position',[.35,.05,0.3,.4],...
    'String','Play','FontSize',12,...
    'Callback',{@button_play,guifig,referencefig,slider_time});
showmaxtemp = uicontrol('Parent',mainpanel,'Style','pushbutton',...
    'Units','normalized','Position',[.025,.1,0.15,.3],...
    'String','MAXTEMP','FontSize',11,...
    'Callback',{@button_maxtemp,guifig,referencefig});
exportvidbutton = uicontrol('Parent',mainpanel,'Style','pushbutton',...
    'Units','normalized','Position',[.675,.1,0.125,.3],...
    'String','save video','FontSize',11,...
    'Callback',{@button_exportvid,guifig,referencefig});
exitbutton = uicontrol('Parent',mainpanel,'Style','pushbutton',...
    'Units','normalized','Position',[.825,.1,0.15,.3],...
    'String','exit + save cfg','FontSize',11,...
    'Callback',{@exportFCN,guifig});
% / main panel

% frame panel
framepanel = uipanel('Title','setup.frames','FontSize',12,...
             'BackgroundColor','white','Position',[.05 .45 .2 .5]);
framepanel_export = uipanel('Parent',framepanel,'Title','export.settings','FontSize',11,...
             'BackgroundColor','white','Position',[0 0 1 .66]);
framestep = uicontrol('Parent',framepanel,'Style','edit',...
    'Units','normalized','Position',[0.5,0.7,0.4,0.2],...
    'String','1',...
    'Callback',{@button_framestep,guifig,referencefig});
uicontrol('Parent',framepanel,'Style','text',...
    'Units','normalized','Position',[0.05,0.65,0.4,0.2],'FontSize',10,...
    'String','stepsize','BackgroundColor',[1 1 1]);
framexportmin = uicontrol('Parent',framepanel_export,'Style','edit',...
    'Units','normalized','Position',[0.5,0.5,0.4,0.35],...
    'String','1',...
    'Callback',{@button_exportmin,guifig,referencefig,slider_time});
uicontrol('Parent',framepanel_export,'Style','text',...
    'Units','normalized','Position',[0.05,0.5,0.4,0.35],'FontSize',10,...
    'String','startframe','BackgroundColor',[1 1 1]);
framexportmax = uicontrol('Parent',framepanel_export,'Style','edit',...
    'Units','normalized','Position',[0.5,0.1,0.4,0.35],...
    'String',num2str(lastframe),...
    'Callback',{@button_exportmax,guifig,referencefig});
uicontrol('Parent',framepanel_export,'Style','text',...
    'Units','normalized','Position',[0.05,0.1,0.4,0.35],'FontSize',10,...
    'String','endframe','BackgroundColor',[1 1 1]);
if plotdefprovided
    framexportmin.String = num2str(data.frames.start);
    framexportmax.String = num2str(data.frames.end);
    framestep.String = num2str(data.frames.step);
end
% / frame panel

% ylim panel
ylimpanel = uipanel('Title','setup.ylims.clims','FontSize',12,...
             'BackgroundColor','white','Position',[.265 .45 .25 .5]);
ylimsinitbtn = getappdata(guifig,'currentylims');
ylimminbutton = uicontrol('Parent',ylimpanel,'Style','edit',...
    'Units','normalized','Position',[0.3,0.7,0.3,0.2],...
    'String',num2str(ylimsinitbtn(1)),...
    'Callback',{@button_ylimmin,guifig,referencefig});
uicontrol('Parent',ylimpanel,'Style','text',...
    'Units','normalized','Position',[0.05,0.675,0.175,0.2],'FontSize',11,...
    'String','min','BackgroundColor',[1 1 1]);
ylimmaxbutton = uicontrol('Parent',ylimpanel,'Style','edit',...
    'Units','normalized','Position',[0.3,0.4,0.3,0.2],...
    'String',num2str(ylimsinitbtn(2)),...
    'Callback',{@button_ylimmax,guifig,referencefig});
uicontrol('Parent',ylimpanel,'Style','text',...
    'Units','normalized','Position',[0.05,0.375,0.175,0.2],'FontSize',11,...
    'String','max','BackgroundColor',[1 1 1]);
ylimresetbutton = uicontrol('Parent',ylimpanel,'Style','pushbutton',...
    'Units','normalized','Position',[0.3,0.1,0.3,0.2],...
    'String','Reset',...
    'Callback',{@button_ylimreset,guifig,referencefig,ylimminbutton,ylimmaxbutton});
templimsinitbtn = getappdata(guifig,'currentclims');
tempminbutton = uicontrol('Parent',ylimpanel,'Style','edit',...
    'Units','normalized','Position',[0.65,0.7,0.3,0.2],...
    'String',num2str(templimsinitbtn(1)),...
    'Callback',{@button_tempmin,guifig,referencefig});
tempmaxbutton = uicontrol('Parent',ylimpanel,'Style','edit',...
    'Units','normalized','Position',[0.65,0.4,0.3,0.2],...
    'String',num2str(templimsinitbtn(2)),...
    'Callback',{@button_tempmax,guifig,referencefig});
tempresetbutton = uicontrol('Parent',ylimpanel,'Style','pushbutton',...
    'Units','normalized','Position',[0.65,0.1,0.3,0.2],...
    'String','Reset',...
    'Callback',{@button_tempreset,guifig,referencefig,tempminbutton,tempmaxbutton});
if plotdefprovided
    ylimminbutton.String = num2str(data.ylims(1),'%0.1f');
    ylimmaxbutton.String = num2str(data.ylims(2),'%0.1f');
    tempminbutton.String = num2str(data.templims(1),'%0.1f');
    tempmaxbutton.String = num2str(data.templims(2),'%0.1f');
end
% / ylim panel

% fig customization stretch panel
customizationpanel = uipanel('Title','fig.customize','FontSize',12,...
             'BackgroundColor','white','Position',[.5325 .45 .2 .5]);
fontplus = uicontrol('Parent',customizationpanel,'Style','pushbutton',...
    'Units','normalized','Position',[0.05,0.75,0.45,0.175],...
    'String','Font +',...
    'Callback',{@button_fontplus,guifig,referencefig});
fontminus = uicontrol('Parent',customizationpanel,'Style','pushbutton',...
    'Units','normalized','Position',[0.5,0.75,0.45,0.175],...
    'String','Font -',...
    'Callback',{@button_fontminus,guifig,referencefig});
ticksplus = uicontrol('Parent',customizationpanel,'Style','pushbutton',...
    'Units','normalized','Position',[0.05,0.5,0.45,0.175],...
    'String','Ticks +',...
    'Callback',{@button_ticksplus,guifig,referencefig});
ticksminus = uicontrol('Parent',customizationpanel,'Style','pushbutton',...
    'Units','normalized','Position',[0.5,0.5,0.45,0.175],...
    'String','Ticks -',...
    'Callback',{@button_ticksminus,guifig,referencefig});
levelsplus = uicontrol('Parent',customizationpanel,'Style','pushbutton',...
    'Units','normalized','Position',[0.05,0.25,0.45,0.175],...
    'String','Levels +',...
    'Callback',{@button_levelsplus,guifig,referencefig});
levelsminus = uicontrol('Parent',customizationpanel,'Style','pushbutton',...
    'Units','normalized','Position',[0.5,0.25,0.45,0.175],...
    'String','Levels -',...
    'Callback',{@button_levelsminus,guifig,referencefig});
stretchdata = uicontrol('Parent',customizationpanel,'Style','edit',...
    'Units','normalized','Position',[0.5,0.025,0.45,0.175],...
    'String','1','Callback',{@button_stretchplus,guifig,referencefig});
uicontrol('Parent',customizationpanel,'Style','text',...
    'Units','normalized','Position',[0.05,0.015,0.45,0.175],'FontSize',10,...
    'String','dataAR','BackgroundColor',[1 1 1]);
if plotdefprovided
   stretchdata.String = num2str(data.stretch_horizontal,2); 
end
% / fig customization stretch panel

% title selector [.7 .45 .2 .25]
title_selector = uibuttongroup('Parent',guifig,'Visible','off',...
                  'Position',[.75 .7 .2 .25],'BackgroundColor',[1 1 1],...
                  'Title','title.selector','FontSize',12,...
                  'SelectionChangedFcn',{@titleselection,referencefig});
uicontrol(title_selector,'Style','radiobutton',...
                  'String','Title #1','BackgroundColor',[1 1 1],...
                  'Position',[10 28 70 30],...
                  'HandleVisibility','off');
uicontrol(title_selector,'Style','radiobutton',...
                  'String','Title #2','BackgroundColor',[1 1 1],...
                  'Position',[10 4 70 30],...
                  'HandleVisibility','off');
title_selector.Visible = 'on'; % visible after child creation
% / title selector and buttons within

% fig vidsettings panel
vidsettings = uipanel('Title','vid.settings','FontSize',12,...
             'BackgroundColor','white','Position',[.75 .45 .2 .25]);
vidfps = uicontrol('Parent',vidsettings,'Style','popupmenu',...
    'Units','normalized','Position',[0.5,0.6,0.45,0.35],...
    'String',{'60','30','25','20','15','10','5'},...
    'Callback',{@button_setfps,guifig,referencefig});
vidfps.Value = 2; % set initial to 30fps
uicontrol('Parent',vidsettings,'Style','text',...
    'Units','normalized','Position',[0.05,0.6,0.45,0.35],'FontSize',10,...
    'String','fps','BackgroundColor',[1 1 1]);
vidcodec = uicontrol('Parent',vidsettings,'Style','popupmenu',...
    'Units','normalized','Position',[0.5,0.1,0.45,0.35],...
    'String',{'mp4','avi'},'Callback',{@button_setcodec,guifig,referencefig});
uicontrol('Parent',vidsettings,'Style','text',...
    'Units','normalized','Position',[0.05,0.1,0.45,0.35],'FontSize',10,...
    'String','codec','BackgroundColor',[1 1 1]);
if plotdefprovided
   stretchdata.String = num2str(data.stretch_horizontal,'%0.1f');
   vidfps.Value = find(ismember(vidfps.String,num2str(data.export.fps)));
   vidcodec.Value = find(ismember(vidcodec.String,num2str(data.export.codec)));
end
% / fig vidsettings panel

% listener executes upon moving the slider instead of after setting the slider
[~] = addlistener(slider_time, 'Value', 'PostSet',...
    @(hObject,eventData) slider_time_callback(hObject,eventData,referencefig,...
                         playbutton,slider_time.Value));
% / listener

% update 1 frame (relevant if plotdefprovided adjusts stuff
drillheataccResults.plot4update(data)

%% now wait for user inputs before outputting the export config spec
uiwait(guifig) % wait for button_export press / fig close
try
    % upon button_export update current data
    output = guidata(referencefig);
    % do not output complete copy of results and remove some other stuff
    output = rmfield(output,{'img'});
    output = rmfield(output,{'labels'});
    output = rmfield(output,{'effective_pulse'});
    output = rmfield(output,{'titledata'});
    % save figure size which may have been adjusted
    output.figsz = [160, 160, output.fig.Position(3:4)];
    % now remove all fields that contain handles, dont need it
    fn = fieldnames(output);
    for i = 1:length(fn)
        if any(ishandle(output.(fn{i}))) && ~isnumeric(output.(fn{i}))
            % weird..ishandle fails if a double == 1 is inside
            output = rmfield(output,(fn{i}));
        end
    end
    % now close / delete actual figure
    delete(referencefig);
catch
    fprintf('Fig %s closed. If you want to save plotfg close through Save/Exit button in GUI\n',uuid)
end

end

function button_exportvid(~,~,guifig,referencefig)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before video export. Execute again.')
    skip = true;
end
if ~skip && isvalid(referencefig)
    % retrieve data
    data = guidata(referencefig);
    % enable export
    data.export.enable = 1;
    % update plot
    drillheataccResults.plot4update(data)
end
end

function button_setfps(hObject,~,guifig,referencefig)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting fps settings. FPS not updated.')
    skip = true;
end
if ~skip && isvalid(referencefig)
    % retrieve data
    data = guidata(referencefig);
    % adjust fps
    data.export.fps = str2double(hObject.String{hObject.Value});
    setappdata(guifig,'vidfps',num2str(data.export.fps));
    % update data
    guidata(referencefig,data)
end
hObject.Value = find(ismember(hObject.String,getappdata(guifig,'vidfps')));
fprintf([guifig.Name,':',32,'fps set to %s\n'],getappdata(guifig,'vidfps'));
reportvidruntime(guifig); % output updated vid export runtime to console
end

function button_setcodec(hObject,~,guifig,referencefig)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting codec settings. Codec not updated.')
    skip = true;
end
if ~skip && isvalid(referencefig)
    % retrieve data
    data = guidata(referencefig);
    % adjust levels
    data.export.codec = hObject.String{hObject.Value};
    setappdata(guifig,'vidcodec',data.export.codec);
    % update data
    guidata(referencefig,data)
end
hObject.Value = find(ismember(hObject.String,getappdata(guifig,'vidcodec')));
fprintf([guifig.Name,':',32,'Codec set to %s\n'],getappdata(guifig,'vidcodec'));
end

function slider_time_callback(~,~,referencefig,playbutton,sliderval)
if isvalid(referencefig)
    data = guidata(referencefig); % retrieve data
    % set frame
    data.frames.currentframe = round(sliderval,0);
    % update plot
    drillheataccResults.plot4update(data)
    % update label
    playbutton.String = sprintf('Frame %i / %i',data.frames.currentframe,data.frames.lastframe);
    % update data
    guidata(referencefig,data)
end
end

function button_play(hObject,~,guifig,referencefig,slider_time)
if ~hObject.Value
    setappdata(guifig,'playing',0);
end
while hObject.Value && isvalid(referencefig)
    if getappdata(guifig,'forcestop')
        hObject.Value = 0;
        setappdata(guifig,'forcestop',0);
        setappdata(guifig,'playing',0);
        break
    end
    setappdata(guifig,'playing',1);
    % retrieve data
    data = guidata(referencefig);
    % advance frame
    data.frames.currentframe = data.frames.currentframe+data.frames.step;
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    guidata(referencefig,data)
    % if lastframe disable play
    % either way update timeslider and label
    if data.frames.currentframe >= data.frames.lastframe
        hObject.Value = 0;
        slider_time.Value = data.frames.lastframe;
        hObject.String = sprintf('Frame %i / %i',data.frames.lastframe,data.frames.lastframe);
    else
        slider_time.Value = data.frames.currentframe;
        hObject.String = sprintf('Frame %i / %i',data.frames.currentframe,data.frames.lastframe);
    end
end
end

function button_stretchplus(hObject,~,guifig,referencefig)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting frame settings. Framestep not updated.')
    skip = true;
end
stretchval = fixnumber(hObject.String,1);
if ~skip && isvalid(referencefig) && ~isempty(stretchval)
    if stretchval < 0.1; stretchval = 0.1; end
    if stretchval > 10; stretchval = 10; end
    setappdata(guifig,'stretchval',stretchval);
    % retrieve data
    data = guidata(referencefig);
    % adjust levels
    data.stretch_horizontal = stretchval;
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    hObject.String = num2str(stretchval); % update editbutton label
    guidata(referencefig,data)
elseif isempty(stretchval)
    hObject.String = num2str(getappdata(guifig,'stretchval'));
end
fprintf([guifig.Name,':',32,'Horizontal data stretch set %1.1f\n'],getappdata(guifig,'stretchval'));
end

function button_fontplus(hObject,~,guifig,referencefig)
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
end
if isvalid(referencefig)
    % retrieve data
    data = guidata(referencefig);
    % advance frame
    data.fontsz = data.fontsz+hObject.Value;
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    guidata(referencefig,data)    
end
end

function button_fontminus(hObject,~,guifig,referencefig)
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
end
if isvalid(referencefig)
    % retrieve data
    data = guidata(referencefig);
    % advance frame
    data.fontsz = data.fontsz-hObject.Value;
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    guidata(referencefig,data)    
end
end

function button_levelsplus(hObject,~,guifig,referencefig)
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
end
if isvalid(referencefig)
    % retrieve data
    data = guidata(referencefig);
    % adjust levels
    if strcmpi(class(data.plot),data.plottypes{1})
        % contourf
        data.contour.levels = data.contour.levels+hObject.Value;
    elseif strcmpi(class(data.plot),data.plottypes{2})
        % surface
        data.surface.levels = data.surface.levels+1;
    end
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    guidata(referencefig,data)    
end
end

function button_levelsminus(hObject,~,guifig,referencefig)
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
end
if isvalid(referencefig)
    % retrieve data
    data = guidata(referencefig);
    % adjust levels
    if strcmpi(class(data.plot),data.plottypes{1})
        % contourf
        data.contour.levels = data.contour.levels-hObject.Value;
    elseif strcmpi(class(data.plot),data.plottypes{2})
        % surface
        data.surface.levels = data.surface.levels-1;
    end
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    guidata(referencefig,data)    
end
end

function button_ticksplus(hObject,~,guifig,referencefig)
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
end
if isvalid(referencefig)
    % retrieve data
    data = guidata(referencefig);
    % adjust levels
    data.numofxticks = data.numofxticks+hObject.Value;
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    guidata(referencefig,data)    
end
end

function button_ticksminus(hObject,~,guifig,referencefig)
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
end
if isvalid(referencefig)
    % retrieve data
    data = guidata(referencefig);
    % adjust levels
    data.numofxticks = data.numofxticks-hObject.Value;
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    guidata(referencefig,data)    
end
end

function button_ylimmin(hObject,~,guifig,referencefig)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting z-limits. Limits not updated.')
    skip = true;
end
maxlims = getappdata(guifig,'ylimsinit');
currentlims = getappdata(guifig,'currentylims');
fixednumber = fixnumber(hObject.String,1);
if ~skip && isvalid(referencefig) && ~isempty(fixednumber)
    if fixednumber < maxlims(1)
        fixednumber = maxlims(1); % dont exceed max
    elseif fixednumber > currentlims(2)
        fixednumber = currentlims(1);
        warning('Lower limit must be smaller than upper limit.');
    end
    % update limits
    updatedlims = [fixednumber, currentlims(2)];
    % update settings in guifig
    setappdata(guifig,'currentylims',updatedlims);
    % retrieve data
    data = guidata(referencefig);
    % adjust ylims
    data.ylims = updatedlims;
    % update plot
    drillheataccResults.plot4update(data)
    % update editbutton label
    hObject.String = num2str(fixednumber);
    % update data
    guidata(referencefig,data)
elseif isempty(fixednumber)
    hObject.String = num2str(currentlims(1));
end
currentlims = getappdata(guifig,'currentylims'); % update current lims
fprintf([guifig.Name,':',32,'Min depth set to %1.1f µm\n'],currentlims(1));
end

function button_ylimmax(hObject,~,guifig,referencefig)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting z-limits. Limits not updated.')
    skip = true;
end
maxlims = getappdata(guifig,'ylimsinit');
currentlims = getappdata(guifig,'currentylims');
fixednumber = fixnumber(hObject.String,1);
if ~skip && isvalid(referencefig) && ~isempty(fixednumber)
    if fixednumber > maxlims(2)
        fixednumber = maxlims(2); % dont exceed max
    elseif fixednumber < currentlims(1)
        fixednumber = currentlims(2);
        warning('Upper limit must be greater than lower limit.');
    end
    % update limits
    updatedlims = [currentlims(1), fixednumber];
    % update settings in guifig
    setappdata(guifig,'currentylims',updatedlims);
    % retrieve data
    data = guidata(referencefig);
    % adjust ylims
    data.ylims = updatedlims;
    % update plot
    drillheataccResults.plot4update(data)
    % update editbutton label
    hObject.String = num2str(fixednumber);
    % update data
    guidata(referencefig,data)
elseif isempty(fixednumber)
    hObject.String = num2str(currentlims(2));
end
currentlims = getappdata(guifig,'currentylims'); % update current lims
fprintf([guifig.Name,':',32,'Max depth set to %1.1f µm\n'],currentlims(2));
end

function button_ylimreset(hObject,~,guifig,referencefig,ylimminbutton,ylimmaxbutton)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting z-limits. Limits not updated.')
    skip = true;
end
maxlims = getappdata(guifig,'ylimsinit');
if ~skip && isvalid(referencefig) && hObject.Value
    % update settings in guifig
    setappdata(guifig,'currentylims',maxlims);
    % retrieve data
    data = guidata(referencefig);
    % adjust ylims
    data.ylims = maxlims;
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    guidata(referencefig,data)
end
currentlims = getappdata(guifig,'currentylims'); % update current lims
ylimminbutton.String = num2str(currentlims(1));
ylimmaxbutton.String = num2str(currentlims(2));
fprintf([guifig.Name,':',32,'Limits reset to %1.1f / %1.1f µm\n'],currentlims(1),currentlims(2));
end

function button_tempmin(hObject,~,guifig,referencefig)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting z-limits. Limits not updated.')
    skip = true;
end
currentlims = getappdata(guifig,'currentclims');
fixednumber = fixnumber(hObject.String,1);
if ~skip && isvalid(referencefig) && ~isempty(fixednumber)
    if fixednumber < 0
        fixednumber = 0; % deltaT < 0 not possible...
    elseif fixednumber > currentlims(2)
        fixednumber = currentlims(1);
        warning('Lower limit must be smaller than upper limit.');
    end
    % update limits
    updatedlims = [fixednumber, currentlims(2)];
    % update settings in guifig
    setappdata(guifig,'currentclims',updatedlims);
    % retrieve data
    data = guidata(referencefig);
    % adjust clims
    data.templims = updatedlims;
    % update plot
    drillheataccResults.plot4update(data)
    % update editbutton label
    hObject.String = num2str(fixednumber);
    % update data
    guidata(referencefig,data)
elseif isempty(fixednumber)
    hObject.String = num2str(currentlims(1));
end
currentlims = getappdata(guifig,'currentclims'); % update current lims
fprintf([guifig.Name,':',32,'Lower temperature display range set to %1.1f K\n'],currentlims(1));
end

function button_tempmax(hObject,~,guifig,referencefig)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting z-limits. Limits not updated.')
    skip = true;
end
currentlims = getappdata(guifig,'currentclims');
fixednumber = fixnumber(hObject.String,1);
if ~skip && isvalid(referencefig) && ~isempty(fixednumber)
    if fixednumber < currentlims(1)
        fixednumber = currentlims(2);
        warning('Upper limit must be greater than lower limit.');
    end
    % update limits
    updatedlims = [currentlims(1), fixednumber];
    % update settings in guifig
    setappdata(guifig,'currentclims',updatedlims);
    % retrieve data
    data = guidata(referencefig);
    % adjust clims
    data.templims = updatedlims;
    % update plot
    drillheataccResults.plot4update(data)
    % update editbutton label
    hObject.String = num2str(fixednumber);
    % update data
    guidata(referencefig,data)
elseif isempty(fixednumber)
    hObject.String = num2str(currentlims(2));
end
currentlims = getappdata(guifig,'currentclims'); % update current lims
fprintf([guifig.Name,':',32,'Upper temperature display range set to %1.1f K\n'],currentlims(2));
end

function button_tempreset(hObject,~,guifig,referencefig,tempminbutton,tempmaxbutton)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting z-limits. Limits not updated.')
    skip = true;
end
maxtemp = getappdata(guifig,'mat_tcrit');
if ~skip && isvalid(referencefig) && hObject.Value
    % update settings in guifig
    setappdata(guifig,'currentclims',[0,maxtemp]);
    % retrieve data
    data = guidata(referencefig);
    % adjust clims
    data.templims = [0,maxtemp];
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    guidata(referencefig,data)
end
currentlims = getappdata(guifig,'currentclims'); % update current lims
tempminbutton.String = '0';
tempmaxbutton.String = num2str(currentlims(2));
fprintf([guifig.Name,':',32,'Temperature display range reset to %1.1f / %4.1f K\n'],currentlims(1),currentlims(2));
end

function button_maxtemp(hObject,~,guifig,referencefig)
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
end
if isvalid(referencefig)
    % retrieve data
    data = guidata(referencefig);
    % adjust levels
    data.frames.showmax = hObject.Value;
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    data.frames.showmax = 0; % reset
    guidata(referencefig,data)
end
end

function titleselection(hObject,eventdata,referencefig)
if getappdata(hObject.Parent,'playing')
    setappdata(hObject.Parent,'forcestop',1);
    warning('Disable PLAY before adjusting title settings. Title not updated.')
end
% get value
switch lower(eventdata.NewValue.String)
    case 'title #1'
        titleidx = 1;
    case 'title #2'
        titleidx = 2;
end
if isvalid(referencefig)
    % retrieve data
    data = guidata(referencefig);
    % adjust levels
    data.titlestr{3} = titleidx;
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    guidata(referencefig,data)
end
end

function button_framestep(hObject,~,guifig,referencefig)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting frame settings. Framestep not updated.')
    skip = true;
end
fixednumber = fixnumber(hObject.String,0);
if ~skip && isvalid(referencefig) && ~isempty(fixednumber)
    if fixednumber < 1; fixednumber = 1; end
    setappdata(guifig,'framestep',fixednumber);
    % retrieve data
    data = guidata(referencefig);
    % adjust levels
    data.frames.step = fixednumber;
    % update plot
    drillheataccResults.plot4update(data)
    % update data
    hObject.String = num2str(fixednumber); % update editbutton label
    guidata(referencefig,data)
elseif isempty(fixednumber)
    hObject.String = num2str(getappdata(guifig,'framestep'));
end
fprintf([guifig.Name,':',32,'Framestepping set to every %i frames\n'],getappdata(guifig,'framestep'));
reportvidruntime(guifig)
end

function button_exportmin(hObject,~,guifig,referencefig,slider_time)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting z-limits. Limits not updated.')
    skip = true;
end
lastrange = getappdata(guifig,'framerange_export');
frame = fixnumber(hObject.String,0);
if ~skip && isvalid(referencefig) && ~isempty(frame)
    if frame >= lastrange(2)
        frame = lastrange(2)-1; % dont exceed max
    elseif frame < 1
        frame = 1;
    end
    % update limits
    newrange = [frame, lastrange(2)];
    % update settings in guifig
    setappdata(guifig,'framerange_export',newrange);
    % retrieve data
    data = guidata(referencefig);
    % adjust startframe
    data.frames.start = frame;
    % move timeslider to correct position
    slider_time.Value = frame;
    data.frames.currentframe = frame;
    % update editbutton label
    hObject.String = num2str(frame);
    % update data
    guidata(referencefig,data)
elseif isempty(frame)
    hObject.String = num2str(lastrange(1));
end
currentlims = getappdata(guifig,'framerange_export'); % update current frames
fprintf([guifig.Name,':',32,'Export set to frame %i through %i\n'],currentlims(1),currentlims(2));
reportvidruntime(guifig); % output updated vid export runtime to console
end

function button_exportmax(hObject,~,guifig,referencefig)
skip = false;
if getappdata(guifig,'playing')
    setappdata(guifig,'forcestop',1);
    warning('Disable PLAY before adjusting z-limits. Limits not updated.')
    skip = true;
end
lastrange = getappdata(guifig,'framerange_export');
lastframe = getappdata(guifig,'lastframe');
frame = fixnumber(hObject.String,0);
if ~skip && isvalid(referencefig) && ~isempty(frame)
    if frame >= lastframe
        frame = lastframe; % dont exceed max
    elseif frame <= lastrange(1)
        frame = lastrange(1)+1; % dont be lower than min
    end
    % update limits
    newrange = [lastrange(1), frame];
    % update settings in guifig
    setappdata(guifig,'framerange_export',newrange);
    % retrieve data
    data = guidata(referencefig);
    % adjust startframe
    data.frames.end = frame;
    % update editbutton label
    hObject.String = num2str(frame);
    % update data
    guidata(referencefig,data)
elseif isempty(frame)
    hObject.String = num2str(lastrange(2));
end
currentlims = getappdata(guifig,'framerange_export'); % update current frames
fprintf([guifig.Name,':',32,'Export set to frame %i through %i\n'],currentlims(1),currentlims(2));
reportvidruntime(guifig); % output updated vid export runtime to console
end

function exportFCN(hObject,~,guifig)
if hObject.Value
    delete(guifig);
end
end

function output = fixnumber(inputstr,decimals)
% takes input string and returns positive number with [decimals]
% if output isempty then bad input -> ignore
% Exclude characters, which are accepted by sscanf:
inputstr(ismember(inputstr,'-+eEgG')) = '';
% Convert to one number and back to a string:
output = abs(round(sscanf(inputstr,'%g',1),decimals));
if isempty(output) || isnan(output) || (output == inf)
    output = [];
    warning('nice try...but inf/nan not allowed.');
end
end

function reportvidruntime(guifig)
% report video metrics to commandline
currentlims = getappdata(guifig,'framerange_export');
fps = str2double(getappdata(guifig,'vidfps'));
fprintf([guifig.Name,':',32,'At %i fps frame %i through %i (step %i) has a playtime of %2.1f seconds.\n'],...
    fps,currentlims(1),currentlims(2),getappdata(guifig,'framestep'),...
    (length(currentlims(1):getappdata(guifig,'framestep'):currentlims(2)))/fps);
end