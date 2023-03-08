function ObjVecTrackingViewer()

% Load the processedDVT to track the run at a certain time window

%%  Initialize GUI creation variables

xAxisRange = 640;
yAxisRange = 480;
buttonWidth = 80;
buttonPadding = 10;
buttonHeight = 25;
overallWindowXSize = 900;
overallWindowYSize = 750;
axesXLoc = 50;
axesYLoc = 200;
buttonRow1Height = 100;
buttonRow2Height = buttonPadding;
buttonColumnXLoc = 715; %Equal padding both sides.
buttonColumnYLoc = axesYLoc;
sliderWidth = 640;

%Use the color generator thing for this.
pen = repmat(distinguishable_colors(10,[0,0,0;1,1,1]),2,1);

%  Create and then hide the GUI as it is being constructed.
overallWindow = figure('Visible','off','Position',...
    [100,100,overallWindowXSize,overallWindowYSize]);

%  Construct the components
hAbort = uicontrol('Style','pushbutton','Callback',@abort_Callback,...
    'String','Reset the GUI!','BackgroundColor','red','Position',...
    [buttonColumnXLoc,buttonColumnYLoc+13*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);

hTimeSlider = uicontrol('Style', 'slider','Min',0,'Max',1,'Value',0,...
    'Position',[axesXLoc,buttonRow2Height,sliderWidth,buttonHeight],...
    'Callback',@timeSlider_Callback);
hWindowTimeSizeLabel = uicontrol('Style','text','String','Time Window Size To Display',...
    'Position',[buttonColumnXLoc,buttonRow2Height+buttonHeight,2*buttonWidth,buttonHeight/2]);
hWindowTimeSize = uicontrol('Style','edit','String',1000,...
    'Position',[buttonColumnXLoc,buttonRow2Height,2*buttonWidth,buttonHeight],...
    'Callback',@windowTimeSize_Callback);

hAxes = axes('Units','Pixels','Position',[axesXLoc,axesYLoc,xAxisRange,yAxisRange]);

%%  GUI Initialization tasks

% Prompt for struct file
[filename, pathname] = uigetfile('*.mat', 'Choose the processed dvt or scored file. (*.mat)');
dvtFileName = fullfile(pathname, filename);
% !!!!!! NOTE !!!!!!! Put in error catch in case they hit cancel here.
% Load dvt file contents.
indRecStruct = load(dvtFileName).indRecStruct;
pixelDvt = indRecStruct.world.processedDVT;

% Sets the mouse back to "ignore-clicks-on-the-grid" mode.
%stateResetter;
%buttonResetter;

%STUB LOAD MOVIE HERE
% [filename, pathname] = uigetfile('*.avi', 'Choose the video file.');
% handles.aviFileName = fullfile(pathname, filename);
% !!!!!! NOTE !!!!!!! Put in error catch in case the hit cancel here.

% Change units to normalized so components resize automatically.
set([hTimeSlider, hWindowTimeSizeLabel, hWindowTimeSize,...
    hAxes,hAbort],'Units','normalized');

% Assign the GUI a name to appear in the window title.
set(overallWindow,'Name','Object Vector Run Tracker');
% Move the GUI to the center of the screen.
movegui(overallWindow,'center');
%Draw first set of data using default values.
drawPoints;
% Make the GUI visible.
guiReady = true;
set(overallWindow,'Visible','on');

%% Callbacks

function abort_Callback(source,eventdata) %#ok<*INUSD>
    % Abort whatever you are doing and return to the normal program
    % state.
    if get(hSelectAndRelabelRun,'Value') && ~isempty(eventBuffer)
        % Write the event you are in the process of relabeling back to
        % spiral events before deleting.
        nEventsToAdd = size(eventBuffer,1);             %#ok<*NODEF>
        eventCount = eventCount + nEventsToAdd;
        if nEventsToAdd ==3
            % Add the run markers to the spiral_events list.
            events(eventCount-2:eventCount,:) = eventBuffer;
            % Paint the tracking file appropriately.
            pixelDvt(events(end-2,1):events(end-1,1),end) = ...
                events(end-2,3);
        else
            % Add the run markers to the spiral_events list.
            events(eventCount-1:eventCount,:) = eventBuffer;
            % Paint the tracking file appropriately.
            pixelDvt(events(end-1,1):events(end,1),end) = ...
                events(end-1,3);
        end
    end

    %stateResetter;
    %buttonResetter;

    set(hCleanRun,'Tag','0');

    set(hMarkRun,'Value',0);
    set(hCleanRun,'Value',0);
    set(hDirtyRun,'Value',0);
    set(hIndecisionPoint,'Value',0);     
    set(hSelectAndRelabelRun,'Value',0);
    set(hSelectAndEraseRun,'Value',0);
    set(hSelectAndEraseEvent,'Value',0);
    set(hPlotOnePathOnly,'Value',0);


    set(hMarkRun,'Enable','on');
    set(hCleanRun,'Enable','off');
    set(hDirtyRun,'Enable','off');
    set(hIndecisionPoint,'Enable','on');

    set(hSelectAndRelabelRun,'Enable','on');
    set(hEraseLast,'Enable','on');
    set(hSelectAndEraseRun,'Enable','on');
    set(hSelectAndEraseEvent,'Enable','on');

    %disableRunButtons;
    set(hNoReward,'Enable','off');
    set(hPlotOnePathOnly,'Enable','on');

    drawPoints;
end

function timeSlider_Callback(source,eventdata)
    % Replot data with new chosen time window.
    drawPoints;
end

function windowTimeSize_Callback(source,eventdata)
    % Replot data with new chosen time window size.
    drawPoints;
end

%% Helper functions

function [firstTimePoint, lastTimePoint, xColumn, yColumn, penColor] = statusSet()
    % Returns the time window being viewed, the columns in the dvt file 
    % of the light being used, and the appropriate pen color for the
    % light.
    currentStepSize = str2double(get(hWindowTimeSize,'String'));
    totalTime = size(pixelDvt,1);
    stepFraction = currentStepSize./totalTime;
    set(hTimeSlider,'SliderStep',[stepFraction,stepFraction*10]);
    firstTimePoint = min(round(get(hTimeSlider,'Value').*totalTime)+1,totalTime);
    lastTimePoint = min(max(firstTimePoint-1+currentStepSize,1),totalTime);

    xColumn = 9;
    yColumn = 10;
    penColor = [0,0.2,0];
end

    function drawPoints()
        %Helper function that does all of the drawing work in the grid.
        
        % Get info about the status of the data in the interface.
        [firstPoint, lastPoint, xLightCol, yLightCol, pen99] = statusSet();
        
        % Additional state info.
        nClicks = str2double(get(hAxes,'Tag'));
        axes(hAxes);
        pen(99,:) = pen99;
        changeStarts = [0;find(diff(pixelDvt(firstPoint:lastPoint,end)))];
        changeEnds = [changeStarts(2:end);lastPoint-firstPoint-1];
        drawnow;
        hold off;
        
        %Actual plotting starts here.
        cla;
        hold on;
        % Some plot formatting.
        grid on;
        set(hAxes,'TickLength',[0 0],'XTickLabel','','YTickLabel','');
        axis([0,640,0,480]);

        % Plots the labeled portions in the correct colors.
        for iSegment = 1:length(changeStarts)
            plot(pixelDvt(firstPoint+changeStarts(iSegment):firstPoint+changeEnds(iSegment),xLightCol),...
                pixelDvt(firstPoint+changeStarts(iSegment):firstPoint+changeEnds(iSegment),yLightCol),'.','Color','black');
        end
        
        % Plots the start and end points for the time segment.
        plot(pixelDvt(firstPoint,xLightCol),...
            pixelDvt(firstPoint,yLightCol),'g+','linewidth',3);
        plot(pixelDvt(lastPoint,xLightCol),...
            pixelDvt(lastPoint,yLightCol),'r+','linewidth',3);
        set(hAxes,'Tag',num2str(nClicks)); %Deals with weird fact that the plot fn resets axes 'Tag' value to ''

    end
end