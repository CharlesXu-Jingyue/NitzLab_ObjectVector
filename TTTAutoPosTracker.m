function TTTAutoPosTracker()
%TTAutoPosTracker  semiautomated Triple T maze behavioral scoring.
%   TTAUTOPOSTRACKER prompts the user to select a processed dvt 
%   struct file and writes a .mat file with a name starting with the name
%   of the dvt file. The user has to do two things in the interface:
%   Mark the beginning and end of each time session.
%   Mark run start/end points - done as "gates" the animal has to run
%   through. Has to be done in order.
% 	The program will then label all of the runs and select them as clean or
% 	dirty based on velocity.
%
%   The program then saves the events, timestamps, pixelDVT and settings 
%   in the struct.
%
%
%   Non-built-in functions called:
%        distinguishable_colors
%
%   Written by Jake Olson, October 2015
%   Last updated by Jake Olson, June 2016 - Changed input/output to structs.


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
hLightChooser = uibuttongroup('Units','pixels','SelectionChangeFcn',@lightChooser_Callback,...
    'Position',[buttonColumnXLoc,buttonRow1Height-buttonHeight,2*buttonWidth,buttonHeight*3]);
hLight1 = uicontrol('Style','Radio','String','Light 1',...
    'Position',[0,3*buttonHeight,2*buttonWidth,buttonHeight],...
    'parent',hLightChooser);
hLight2 = uicontrol('Style','Radio','String','Light 2',...
    'Position',[0,2*buttonHeight,2*buttonWidth,buttonHeight],...
    'parent',hLightChooser);
hLight3 = uicontrol('Style','Radio','String','Light 3 (Avg)',...
    'Position',[0,buttonHeight,2*buttonWidth,buttonHeight],...
    'parent',hLightChooser);
hLight4 = uicontrol('Style','Radio','String','Light 4 (Mashup)',...
    'Position',[0,0,2*buttonWidth,buttonHeight],...
    'parent',hLightChooser);

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

% Initialize variables - many are used in local and nested functions.
% Prompt for dvt and avi files. - Replaced when TTTrackingPreprocessor
[filename, pathname] = uigetfile('*.mat', 'Choose the processed dvt file. (*.mat)');
dvtFileName = fullfile(pathname, filename);
% !!!!!! NOTE !!!!!!! Put in error catch in case they hit cancel here.
% Load dvt file contents.
load(dvtFileName, 'indRecStruct');

% Initialize variables so they exist across functions.
events = zeros(0,3);
pixelDvt = indRecStruct.processedDVT;
pixelDvt(:,end+1) = 99; %The last row will be updated with run label values.
timeStamps = [];
objective = 0;
watchKeysFlag = 0;
clickLocation = [];
dvtRow = 0;
nClicks = 0;
instructions = 'Mark the beginning of the platform.'; % 1st instruction.
pathZones = []; %minx miny maxx maxy
dvtZoneMarked = [];
autoRuns = [];

% Save original pointer type to be used in mouse hovering
originalPointerType = get(gcf,'Pointer');
% Set function to run on close.
set(overallWindow, 'DeleteFcn', @hCloseFunction);
set(overallWindow,'WindowKeyPressFcn',@keyPressCallback);
% Set event list print format.
format short g;

% Change units to normalized so components resize automatically.
set([overallWindow,...
    hLightChooser, hLight1, hLight2, hLight3,hLight4,...
    hTimeSlider, hWindowTimeSizeLabel, hWindowTimeSize,...
    hAxes], 'Units','normalized');

% Assign the GUI a name to appear in the window title.
set(overallWindow,'Name','TT Position Tracker');
% Move the GUI to the center of the screen.
movegui(overallWindow,'center');
%Draw first set of data using default values.
drawPoints;
% Make the GUI visible.
guiReady = true;
set(overallWindow,'Visible','on');
set(gcf,'WindowButtonMotionFcn',@detectMouseHover);

%% Three main Objectives
% Mark Time Sessions
% objective = 0;
% Mark run start/end points
% objective = 1;
% Label individual runs
% runFinder - no input needed.

%%  Callbacks for autoPosTracker2015
% Callbacks automatically have access to variables nested at a lower level.
    function lightChooser_Callback(source,eventdata)
        % Replot data with new chosen light.
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

%% Mouse hover and click callbacks
    function detectMouseHover(handles,object, eventdata)
        % Handles the appearance of the mouse and what happens when a
        % click is made.
        if ~guiReady
            return;
        end
        
        newAxesPos = estimateNewPosition;
        currentLocInAxes = get(hAxes,'CurrentPoint');
        cursorInAxes = currentLocInAxes(1,1)>= 1 && currentLocInAxes(1,1)<=xAxisRange...
            && currentLocInAxes(1,2)>= 1  && currentLocInAxes(1,2) <= yAxisRange;
        
        if cursorInAxes
            set(gcf,'Pointer','cross');
            set(gcf,'WindowButtonDownFcn',@detectMousePressRecord);
        else
            set(gcf,'Pointer',originalPointerType);
            % If doing nothing, you need to tell it to do nothing, or else
            % the other condition for 'WindowButtonDownFcn' will be exceuted
            %  from the previous mouse position and that is not desirable.
            set(gcf,'WindowButtonDownFcn','');
        end
    end

    function detectMousePressRecord(handles,object, eventdata)
        % Handles mouse clicks - this is the function that does all of the
        % work - any time a click is recorded the state of the buttons must
        % be determined and the appropriate action taken.
        
        % Get the location of the mouse at time of the click.
        clickLocation = get(hAxes,'CurrentPoint');
        clickLocation = clickLocation(1,1:2);
        % Get info about the status of the data in the interface.
        [firstPoint, lastPoint, xLightCol, yLightCol] = statusSet();
        
        nClicks = nClicks+1;
        
        % Determine plotted points.
        allLinesPlotted = [firstPoint:lastPoint]';
        
        % Determine point currently plotted that click is closest to.
        runDist = dist(clickLocation(1,1:2),pixelDvt(allLinesPlotted,[xLightCol,yLightCol])');
        [~, ind] = min(runDist);
        % this is the row index of the user-selected point, with
        % respect to the full 'pos' matrix.
        dvtRow = allLinesPlotted(ind(1));
        
        % Process click. Three general cases - selecting a run/marking an
        % event, labeling starts and end regions, and marking session
        % times.
        
        switch objective
            case 0
                sessionMarker;
            case 1
                zoneMarker;
            case 2
                runMarker;
        end
    end

    function keyPressCallback(source,eventdata)
        % determine the key that was pressed
%         if strcmp(eventdata.Key,'t')
%             % Save timestamps.
%             save('autoPosTimeStamps', 'timeStamps');
%             close all;
%             return;
%         end
        if watchKeysFlag
            keyPressed = eventdata.Key;
            if objective == 0 & nClicks == 4 & strcmpi(keyPressed,'return')
                % no track rotation condition
                objective = 1;
                nClicks = 0;
                instructions = 'PATHS 1-8 - Mark the start line with two clicks.'; %091218
                moveTimes(1);
                watchKeysFlag = 0;
            end
            drawPoints;
        end
    end

%% GUI close callback
function hCloseFunction(handles,object,eventdata)
%         % Save spiral_events labeled and the timestamps.
%         save('autoPosTrackerTMPSAVEFILE', 'events','timeStamps');
 end

%%  Program Control for autoPosTracker
    function sessionMarker
        timeStamps(nClicks) = dvtRow;
        switch nClicks
            case 0
                instructions = 'Mark the beginning of the platform.';
            case 1
                instructions = 'Mark the end of the platform.';
            case 2
                instructions = 'Mark the beginning of the 1st session.';
            case 3
                instructions = 'Mark the end of the 1st session.';
            case 4
                instructions = 'Mark the beginning of the 2nd session.';
                watchKeysFlag = 1;
            case 5
                instructions = 'Mark the end of the 2nd session.';
                watchKeysFlag = 0;
        end
        
        if nClicks == 6
            objective = 1;
            nClicks = 0;
            instructions = 'PATHS 1-8 - Mark the start line with two clicks.'; %091218
            moveTimes(1);
        end
        drawPoints;
    end

    function zoneMarker
        switch nClicks
            case {0,1}
                instructions = 'PATHS 1-8 - Mark the start line with two clicks.';
            case {2,3}
                instructions = 'PATH 1 - Mark the end line with two clicks.';
            case {4,5}
                instructions = 'PATH 2 - Mark the end line with two clicks.';
            case {6,7}
                instructions = 'PATH 3 - Mark the end line with two clicks.';
            case {8,9}
                instructions = 'PATH 4 - Mark the end line with two clicks.';
            case {10,11}
                instructions = 'PATH 5 - Mark the end line with two clicks.';
            case {12,13}
                instructions = 'PATH 6 - Mark the end line with two clicks.';
            case {14,15}
                instructions = 'PATH 7 - Mark the end line with two clicks.';
            case {16,17}
                instructions = 'PATH 8 - Mark the end line with two clicks.';
            case {18,19}
                instructions = 'PATH 9 - Mark the start line with two clicks.';
            case {20,21}
                instructions = 'PATH 9 - Mark the end line with two clicks.';
            case {22,23}
                instructions = 'PATH 10 - Mark the start line with two clicks.';
            case {24,25}
                instructions = 'PATH 10 - Mark the end line with two clicks.';
            case {26,27}
                if length(timeStamps) >10 %4 t o10 091218 (10 paths?) 
                    moveTimes(2);
                end
                instructions = 'PATHS 11-14 - Mark the start line with two clicks.';
            case {28,29}
                instructions = 'PATH 11 - Mark the end line with two clicks.';
            case {30,31}
                instructions = 'PATH 12 - Mark the end line with two clicks.';
            case {32,33}
                instructions = 'PATH 13 - Mark the end line with two clicks.';
            case {34,35}
                instructions = 'PATH 14 - Mark the end line with two clicks.';
            case {36,37}
                instructions = 'PATH 19 - Mark the start line with two clicks.';
            case {38,39}
                instructions = 'PATH 19 - Mark the end line with two clicks.';
            case {40,41}
                instructions = 'PATH 20 - Mark the start line with two clicks.';
            case {42,43}
                instructions = 'PATH 20 - Mark the end line with two clicks.';
            otherwise
        end
        pathZones(nClicks,:) = clickLocation;
        
        if (length(timeStamps) == 8 && nClicks == 27) || nClicks == 27 %091218 changed from 4 to 8 nclicks now only for session 1
            runFinder;
            save(strcat(dvtFileName(1:end-16),'PostAutoEvents'),'indRecStruct');
            close(overallWindow);
        else
            drawPoints;
        end
    end

    function runFinder
        % Finds all of the runs that we are going to label.
        
        % Calculate the cross points
        for iPathZone = 1:floor(size(pathZones,1)/2)
            minX = min(pathZones(2*iPathZone-1:2*iPathZone,1));
            minY = min(pathZones(2*iPathZone-1:2*iPathZone,2));
            maxX = max(pathZones(2*iPathZone-1:2*iPathZone,1));
            maxY = max(pathZones(2*iPathZone-1:2*iPathZone,2));
            
            distX = maxX-minX;
            distY = maxY-minY;
            
            % Get light columns
            [~,~, xLightCol, yLightCol,~] = statusSet();
            
            % Label all gate crossings
            if distX > distY %Horizontal line.
                lineVal = (minY+maxY)/2; % Y value of line
                dvtZoneMarked(:,iPathZone) = ...
                    pixelDvt(:,xLightCol) >= minX &...
                    pixelDvt(:,xLightCol) <= maxX &...
                    ((pixelDvt(:,yLightCol) >= lineVal &...
                    [pixelDvt(2:end,yLightCol);1] <= lineVal) |...
                    (pixelDvt(:,yLightCol) <= lineVal &...
                    [pixelDvt(2:end,yLightCol);1] >= lineVal));
            else %Vertical line
                lineVal = (minX+maxX)/2; % X value of line
                dvtZoneMarked(:,iPathZone) = ...
                    pixelDvt(:,yLightCol) >= minY &...
                    pixelDvt(:,yLightCol) <= maxY &...
                    ((pixelDvt(:,xLightCol) >= lineVal &...
                    [pixelDvt(2:end,xLightCol);1] <= lineVal) |...
                    (pixelDvt(:,xLightCol) <= lineVal &...
                    [pixelDvt(2:end,xLightCol);1] >= lineVal));
            end
        end
        %Erase any gate crossings from outside the appropriate time
        %sessions.
        if length(timeStamps) == 10 %091218
            sessFilter = zeros(length(pixelDvt),26); %091218 changed from 18 to 26 (max click number)? 
            sessFilter(timeStamps(3):timeStamps(4),1:13) = 1; %091218 changed 1:9 to 1:13 number of gates
            sessFilter(timeStamps(5):timeStamps(6),14:26) = 1; %091218 changed 10:18 to 18:26
            dvtZoneMarked = dvtZoneMarked & sessFilter;
        end
        
        
        % label potential runs.
        for iSample = timeStamps(3):timeStamps(4)
            if dvtZoneMarked(iSample,1) % 1st session runs
                [runType, line] = find(dvtZoneMarked(iSample+1:timeStamps(4),:)',1);
                if ~isempty(runType)
                    if any(runType == [2:9]) %091218 changed 2:5 to 2:9 +4 for 4 new paths
                        autoRuns = [autoRuns; iSample,line+iSample,runType-1];
                    end
                end
            elseif dvtZoneMarked(iSample,10) %091218 6 to 10 gate? 
                [runType, line] = find(dvtZoneMarked(iSample+1:timeStamps(4),:)',1);
                if ~isempty(runType)
                    if runType == 11 %091218 changed to fit pattern 7 to 11
                        autoRuns = [autoRuns; iSample,line+iSample,9];
                    end
                end
            elseif dvtZoneMarked(iSample,12) %091218 8 to 12 gate? 
                [runType, line] = find(dvtZoneMarked(iSample+1:timeStamps(4),:)',1);
                if ~isempty(runType)
                    if runType == 13 %091218 changed to fit pattern 9 to 13 
                        autoRuns = [autoRuns; iSample,line+iSample,10];
                    end
                end
            end
        end
        if length(timeStamps) == 10;
            for iSample = timeStamps(5):timeStamps(6)
                if dvtZoneMarked(iSample,10) % 1st session runs
                    [runType, line] = find(dvtZoneMarked(iSample+1:timeStamps(6),:)',1);
                    if ~isempty(runType)
                        if any(runType == [19:22]) %091218 changed with click #s for route 9 11:14 to 19:22
                            autoRuns = [autoRuns; iSample,line+iSample,runType];
                        end
                    end
                elseif dvtZoneMarked(iSample,15)
                    [runType, line] = find(dvtZoneMarked(iSample+1:timeStamps(6),:)',1);
                    if ~isempty(runType)
                        if runType == 24 %091218 kept with pattern 16 to 24
                            autoRuns = [autoRuns; iSample,line+iSample,19];
                        end
                    end
                elseif dvtZoneMarked(iSample,17)
                    [runType, line] = find(dvtZoneMarked(iSample+1:timeStamps(6),:)',1);
                    if ~isempty(runType)
                        if runType == 26  %091218 kept with pattern 18 to 26
                            autoRuns = [autoRuns; iSample,line+iSample,20];
                        end
                    end
                end
            end
        end
        
        %Label clean/dirty and convert to spiral events.
        goodLight = find(min(sum(indRecStruct.samplesUnfilled)));
        minSpeedAllowedCleanRun = 10; % in pixels/sec becasue that is what the velSmoothed is in.
        for iAutoRun = 1:length(autoRuns)
            if any(indRecStruct.velSmoothed(autoRuns(iAutoRun,1):autoRuns(iAutoRun,2),1,goodLight)<minSpeedAllowedCleanRun);
                events = [events;...
                    [pixelDvt(autoRuns(iAutoRun,1),1:2),-1*autoRuns(iAutoRun,3)];...
                    [pixelDvt(autoRuns(iAutoRun,2),1:2),-1*(autoRuns(iAutoRun,3)+100)]];
            else
                events = [events;...
                    [pixelDvt(autoRuns(iAutoRun,1),1:2),autoRuns(iAutoRun,3)];...
                    [pixelDvt(autoRuns(iAutoRun,2),1:2),autoRuns(iAutoRun,3)+100]];
            end
        end
        
        indRecStruct.sessionTimeStamps = timeStamps;
        indRecStruct.events = events;
        indRecStruct.pixelDVT = pixelDvt;
        indRecStruct.minSpeedAllowed = minSpeedAllowedCleanRun;
    end

    

%%  Utility functions for autoPosTracker
    function newAxesPos = estimateNewPosition()
        %% This function is used to get the new position of the uipanel, to be used for GUIs that are Resizable.
        guiPosition = get(overallWindow,'Position');
        axesPosition = get(hAxes,'Position');
        
        newAxesPos(1) = guiPosition(3)*axesPosition(1);
        newAxesPos(2) = guiPosition(4)*axesPosition(2);
        newAxesPos(3) = guiPosition(3)*axesPosition(3);
        newAxesPos(4) = guiPosition(4)*axesPosition(4);
    end

    function moveTimes(session)
        % Changes the times shown to the appropriate window.
        set(hTimeSlider,'Value',timeStamps(session*2+1)/length(pixelDvt));
        set(hWindowTimeSize,'String',num2str(timeStamps(session*2+2)-timeStamps(session*2+1)));
        drawPoints;
    end

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
        
        currentLight = get(hLightChooser,'SelectedObject');
        if currentLight == hLight1
            xColumn = 3;
            yColumn = 4;
            penColor = [0.4,0.4,0.4];
        elseif currentLight == hLight2
            xColumn = 5;
            yColumn = 6;
            penColor = [0.2,0.2,0.2];
        elseif currentLight == hLight3
            xColumn = 7;
            yColumn = 8;
            penColor = [0,0,0];
        else
            xColumn = 9;
            yColumn = 10;
            penColor = [0,0.2,0];
        end
    end

    function drawPoints
        % Helper function that does all of the drawing work in the grid.
        
        % Get info about the status of the data in the interface.
        [firstPoint, lastPoint, xLightCol, yLightCol, pen99] = statusSet();
        
        % Additional state info.
        axes(hAxes);
        pen(99,:) = pen99;
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
        plot(pixelDvt(firstPoint:lastPoint,xLightCol),...
            pixelDvt(firstPoint:lastPoint,yLightCol),'.',...
            'Color',pen(abs(pixelDvt(lastPoint,end)),:));
        
        % Plots the start and end points for the time segment.
        plot(pixelDvt(firstPoint,xLightCol),...
            pixelDvt(firstPoint,yLightCol),'g+','linewidth',3);
        plot(pixelDvt(lastPoint,xLightCol),...
            pixelDvt(lastPoint,yLightCol),'r+','linewidth',3);
        
        % Plot timeStamp points
        plottedTimeStamps = find(timeStamps >= firstPoint & timeStamps <= lastPoint);
        if plottedTimeStamps
            for iStamp = 1:length(plottedTimeStamps)
                text(pixelDvt(timeStamps(plottedTimeStamps(iStamp)),xLightCol),...
                    pixelDvt(timeStamps(plottedTimeStamps(iStamp)),yLightCol),...
                    'O','FontSize',15,'FontWeight','bold','Color','b');
            end
        end
        
        % Plot run start/stop lines.
        if objective == 1 && nClicks < 26 %091218 changed with max clicks
            if size(pathZones,1) > 1
                for iPathZone = 1:floor(size(pathZones,1)/2)
                    plot(pathZones(2*iPathZone-1:2*iPathZone,1)',...
                        pathZones(2*iPathZone-1:2*iPathZone,2)','cd-');
                end
            end
            if mod(nClicks,2)
                plot(pathZones(end,1),pathZones(end,2),'gd')
            end
        elseif objective == 1
            secondSessionOffset = 26; %091218 changed with max clicks
            if size(pathZones,1) > 1 + secondSessionOffset
                for iPathZone = 1:floor((size(pathZones,1)-secondSessionOffset)/2)
                    plot(pathZones(2*iPathZone-1+secondSessionOffset:2*iPathZone+secondSessionOffset,1)',...
                        pathZones(2*iPathZone-1+secondSessionOffset:2*iPathZone+secondSessionOffset,2)','cd-');
                end
            end
            if mod(nClicks,2)
                plot(pathZones(end,1),pathZones(end,2),'gd')
            end
        end
        % Write instructions.
        text(0,yAxisRange+10,instructions,'FontSize',15,'FontWeight','bold','Color','k');
    end

end
