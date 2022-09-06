function plot4update(data)
% disable warning
warning('off','MATLAB:contour:ConstantData')
% check if we need to change number of xticks
if (length(data.ax.XTick)+1) ~= (data.numofxticks*2)
    updatexticks()
end
% check if we need to update ylimits
if any(data.ax.YLim ~= data.ylims)
    ylim(data.ax,data.ylims)
end
% check if we need to update colorbar limits
if any(data.cbar.Limits ~= data.templims)
    caxis(data.ax,data.templims);
    % in that case we alos need to update cmap if contourf
    if strcmpi(class(data.plot),data.plottypes{1})
        % if its contourf
        updatecontourfcmap();
    elseif strcmpi(class(data.plot),data.plottypes{2})
        % else its surface
        updatecbarticks();
    end
end
% check if we need to change colormap / contourlevels (pressing ticks when contourf)
if strcmpi(class(data.plot),data.plottypes{1}) % check if contour
    if length(data.plot.LevelList) ~= data.contour.levels
        updatecontourfcmap();
    end
end
% check if we need to update surface levels (pressing ticks when surfplot)
if strcmpi(class(data.plot),data.plottypes{2})
    if length(data.cbar.Ticks) ~= data.surface.levels
        updatecbarticks();
    end
end
% check if user request video output
if data.export.enable
    time_export = tic;
    vidwriter = startrecord();
    range = data.frames.start:data.frames.step:data.frames.end;
    len = length(range);
    fprintf('Exporting frames %i through %i in steps of %i\n',...
        data.frames.start,data.frames.end,data.frames.step);
end
%% actual stuff starts here
done = false; % init
frame = 1; % init
while ~done
    if data.export.enable
        data.frames.currentframe = range(frame);
    else
        done = true; % only one frame
    end
    i = data.frames.actualframe;
    % redraw if actualframe ~= currentframe
    if data.frames.showmax || (data.frames.actualframe ~= data.frames.currentframe)
        if ~data.frames.showmax
            if data.frames.currentframe >= data.frames.lastframe
                i = data.frames.lastframe;
            else
                i = data.frames.currentframe;
            end
            img = data.img.zdata(:,:,i);
            hslen = data.img.hslen(i);
        elseif data.frames.showmax
            img = data.img.zdatamax;
            hslen = data.img.hslenmax;
        end
        % prepare image
        im = [flip(img(:,2:end),2), img(:,:)];
        % update plot
        data.plot.ZData = im;
        % update title
        if ~data.frames.showmax
            titleidx = data.titlestr{3};
            data.title = title(data.ax,sprintf(data.titlestr{titleidx},data.titledata{1}(i),...
                data.titledata{2}(i),data.titledata{3}(i),data.titledata{4}(i),...
                data.titledata{5}(i),data.img.hslen(i)),'FontSize',data.fontsz(1));
        elseif data.frames.showmax
            data.title = title(data.ax,'MaxT');
        end
        % update patch
        if data.isdrilling(1)
            data.patch.YData(3) = hslen;
        end
        % highlight contour line
        drawnow % below REQUIRES drawnow
        if strcmpi(class(data.plot),data.plottypes{1})
            % below higlights the MAX isoline in white IFF its associated
            % temperature exceeds the critical temperature of the material
            for j = 1:length(data.plot.EdgePrims)
                data.plot.EdgePrims(j).ColorData = uint8([0;0;0;0]);
            end
            if ~isempty(data.plot.EdgePrims) && (max(im(:)) > data.contour.isolines(end))
                set(data.plot.EdgePrims(end),'LineStyle','solid','LineWidth',2.5)
                data.plot.EdgePrims(end).ColorData = uint8([255;255;255;255]);
            end
        end
    end
    % update annotation
    actualpos = plotboxpos(data.ax);
    offset = (42+1.325*(data.fontsz(1)-14))*data.ax.Position(4)/data.fig.Position(4);
    data.annotation.Position = [actualpos(1), actualpos(2)+offset, 0.5, 0];
    if ~data.frames.showmax
        data.annotation.String = sprintf('%4.0f',data.effective_pulse(i));
    elseif data.frames.showmax
        data.annotation.String = '';
    end
    data.annotation.FontSize = data.fontsz(1);
    % update fontsizes
    data.cbar.FontSize = data.fontsz(2); data.labels.cbar.FontSize = data.fontsz(1);
    data.ax.XAxis.FontSize = data.fontsz(1); data.labels.x.FontSize = data.fontsz(1)+2;
    data.ax.YAxis.FontSize = data.fontsz(1); data.labels.y.FontSize = data.fontsz(1)+2;
    % update data aspect ratio (horizontal stretch)
    daspect(data.ax,[1 data.stretch_horizontal 1]) % set aspect ratio, allow for stretching (2)
    % force redraw
    drawnow, warning('on','MATLAB:contour:ConstantData')
    % advance frame index when exporting video + write frame
    if data.export.enable
        if (frame+1) < len
            frame = frame+1;
            % progress ocnsole
            if mod(frame,floor(length(range)/7)) == 0
                current_time = toc(time_export);
                time_remain = round((current_time/frame)*len-current_time,0);
                fprintf('est. ~%1.1f min [%1.0f%%] -- ',time_remain/60,100*frame/len);
            end
        else
            done = true; fprintf('\n');
        end
        % write one frame
        writeVideo(vidwriter, getframe(data.fig));
    end
end
% reset stuff if video was exported
if data.export.enable
    close(vidwriter); fprintf('Video saved successfully.\n')
    set(data.fig,'MenuBar','figure'); set(data.fig,'ToolBar','figure'); % enable again
    data.export.enable = 0;
end

% IMPORTANT: update / write guidata
data.frames.actualframe = i;
guidata(data.fig,data);

%% internal functions
    function vidwriter = startrecord()
        % fix image size mod8 for vidcodec
        figsz = (data.fig.Position(3:4)-mod(data.fig.Position(3:4),8));
        data.fig.Position(3:4) = figsz;
        % disable fig bars
        set(data.fig,'MenuBar','none'); set(data.fig,'ToolBar','none');
        outputfolder = fullfile(pwd,'results');
        if ~isfolder(outputfolder), mkdir(outputfolder); end
        fullfilename = fullfile(outputfolder,strcat(data.UUID,'_',generateID(3),'.',data.export.codec));
        fprintf('Output video: %s\n',fullfilename);
        switch data.export.codec
            case 'mp4'
                vidwriter = VideoWriter(fullfilename,'MPEG-4');
                vidwriter.Quality = 100; % always 100%
            case 'avi'
                vidwriter = VideoWriter(fullfilename,'Uncompressed AVI');
        end
        vidwriter.FrameRate = data.export.fps;
        open(vidwriter);
    end
    function updatecbarticks()
        ticks_tmp = round(linspace(data.templims(1),data.templims(2),data.surface.levels),-1);
        if length(unique(ticks_tmp)) ~= data.surface.levels
            % can't have copies of ticks at the same temp -> round(n,0) instead of round(n,-1)
            ticks_tmp = round(linspace(data.templims(1),data.templims(2),data.surface.levels),0);
        end
        data.cbar.Ticks = ticks_tmp;
    end

    function updatecontourfcmap()
        % used to update contour(f) isolines, colorbar ticks, colormap
        data.contour.isolines = round(linspace(data.templims(1),data.templims(2),data.contour.levels),0);
        customcmap = interp2(jet(255),1:3,linspace(1,255,data.contour.levels-1).');
        % update contour(f) isolines
        data.plot.LevelList = data.contour.isolines;
        % colorbar ticks
        data.cbar.Ticks = data.contour.isolines;
        % apply correct colormap
        data.ax.Colormap = customcmap;
    end
    function updatexticks()
        % updates x ticks upon request
        xticksvect = round(linspace(0,data.img.radiusvect(end),data.numofxticks),0);
        xticksvect = [-flip(xticksvect(2:end)), xticksvect];
        data.ax.XTick = xticksvect;
        data.ax.XTickLabel = num2cell(abs(xticksvect.')); % no negative radii
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