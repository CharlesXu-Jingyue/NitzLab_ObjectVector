function TTPosTracker()
%TTPosTracker  Triple T maze behavioral scoring interface.
%   TTPosTracker is a behavioral scoring interface for dvt files.
%   To be used after AutoPosTracker to check and correct the scoring.
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
hAbort = uicontrol('Style','pushbutton','Callback',@abort_Callback,...
    'String','Reset the GUI!','BackgroundColor','red','Position',...
    [buttonColumnXLoc,buttonColumnYLoc+13*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);
hSeeEventList = uicontrol('Style','pushbutton','Callback',@seeEventList_Callback,...
    'String','Print Event List','Position',...
    [buttonColumnXLoc,buttonColumnYLoc+12*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);
hSaveEventList = uicontrol('Style','pushbutton','Callback',@saveEventList_Callback,...
    'String','Save Events','BackgroundColor','Green','Position',...
    [buttonColumnXLoc,buttonColumnYLoc+11*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);

hMarkRun = uicontrol('Style','pushbutton','Callback',@markRun_Callback,...
    'String','Mark Run','Position',...
    [buttonColumnXLoc,buttonColumnYLoc+9*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);
hCleanRun = uicontrol('Style','togglebutton','Callback',@cleanRun_Callback,...
    'String','Clean Run!','Tag','0','Position',...
    [buttonColumnXLoc,buttonColumnYLoc+8*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);
hDirtyRun = uicontrol('Style','togglebutton','Callback',@dirtyRun_Callback,...
    'String','Dirty Run','Position',...
    [buttonColumnXLoc,buttonColumnYLoc+7*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);
hIndecisionPoint = uicontrol('Style','togglebutton','Callback',@indecisionPoint_Callback,...
    'String','Mark Indecision Point','Position',...
    [buttonColumnXLoc,buttonColumnYLoc+6*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);


hSelectAndRelabelRun = uicontrol('Style','togglebutton','Callback',@selectAndRelabelRun_Callback,...
    'String','Select and Relabel Run','Position',...
    [buttonColumnXLoc,buttonColumnYLoc+4*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);


hEraseLast = uicontrol('Style','pushbutton','Callback',@eraseLastEvent_Callback,...
    'String','Erase Last Run/Event','Position',...
    [buttonColumnXLoc,buttonColumnYLoc+2*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);
hSelectAndEraseRun = uicontrol('Style','togglebutton','Callback',@selectAndEraseRun_Callback,...
    'String','Select and Erase Run','Position',...
    [buttonColumnXLoc,buttonColumnYLoc+1*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);
hSelectAndEraseEvent = uicontrol('Style','togglebutton','Callback',@selectAndEraseEvent_Callback,...
     'String','Select and Erase Marker','Position',...
     [buttonColumnXLoc,buttonColumnYLoc+0*(buttonHeight+buttonPadding),buttonWidth*2,buttonHeight]);


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

hAddedRunLabels = uicontrol('Style','togglebutton','Callback',@addedRunLabels_Callback,...
    'String','Change buttons to 11-20','Position',...
    [axesXLoc,...
    buttonRow1Height-(buttonHeight+buttonPadding),...
    buttonWidth*2,buttonHeight]);
hPlotOnePathOnly = uicontrol('Style','togglebutton','Callback',@plotOnePathOnly_Callback,...
    'String','Plot only one path''s runs','Position',...
    [axesXLoc+xAxisRange-2*buttonWidth,...
    buttonRow1Height-(buttonHeight+buttonPadding),...
    buttonWidth*2,buttonHeight]);

hRunChooser = uibuttongroup('Units','pixels','SelectionChangeFcn',@runChooser_Callback,...
    'Position',[axesXLoc,buttonRow1Height,xAxisRange,...
    3*buttonPadding+2*buttonHeight]);
hRun1 = uicontrol('Style','togglebutton','parent',hRunChooser,...
    'String','Run 1','Tag','1','BackgroundColor',pen(1,:),...
    'Position',[20,2*buttonPadding+buttonHeight,buttonWidth,buttonHeight]);
hRun2 = uicontrol('Style','togglebutton','parent',hRunChooser,...
    'String','Run 2','Tag','2','BackgroundColor',pen(2,:),...
    'Position',[120,2*buttonPadding+buttonHeight,buttonWidth,buttonHeight]);
hRun3 = uicontrol('Style','togglebutton','parent',hRunChooser,...
    'String','Run 3','Tag','3','BackgroundColor',pen(3,:),...
    'Position',[220,2*buttonPadding+buttonHeight,buttonWidth,buttonHeight]);
hRun4 = uicontrol('Style','togglebutton','parent',hRunChooser,...
    'String','Run 4','Tag','4','BackgroundColor',pen(4,:),...
    'Position',[320,2*buttonPadding+buttonHeight,buttonWidth,buttonHeight]);
hRun5 = uicontrol('Style','togglebutton','parent',hRunChooser,...
    'String','Run 5','Tag','5','BackgroundColor',pen(5,:),...
    'Position',[420,2*buttonPadding+buttonHeight,buttonWidth,buttonHeight]);
hRun6 = uicontrol('Style','togglebutton','parent',hRunChooser,...
    'String','Run 6','Tag','6','BackgroundColor',pen(6,:),...
    'Position',[20,buttonPadding,buttonWidth,buttonHeight]);
hRun7 = uicontrol('Style','togglebutton','parent',hRunChooser,...
    'String','Run 7','Tag','7','BackgroundColor',pen(7,:),...
    'Position',[120,buttonPadding,buttonWidth,buttonHeight]);
hRun8 = uicontrol('Style','togglebutton','parent',hRunChooser,...
    'String','Run 8','Tag','8','BackgroundColor',pen(8,:),...
    'Position',[220,buttonPadding,buttonWidth,buttonHeight]);
hRun9 = uicontrol('Style','togglebutton','parent',hRunChooser,...
    'String','Run 9','Tag','9','BackgroundColor',pen(9,:),...
    'Position',[320,buttonPadding,buttonWidth,buttonHeight]);
hRun10 = uicontrol('Style','togglebutton','parent',hRunChooser,...
    'String','Run 10','Tag','10','BackgroundColor',pen(10,:),...
    'Position',[420,buttonPadding,buttonWidth,buttonHeight]);
hNoReward = uicontrol('Style','togglebutton','parent',hRunChooser,...
    'String','No Reward Given','Tag','No Reward','BackgroundColor','white',...
    'Position',[520,buttonPadding,buttonWidth,buttonHeight]);


hAxes = axes('Units','Pixels','Position',[axesXLoc,axesYLoc,xAxisRange,yAxisRange]);


%%  GUI Initialization tasks
% Prompt for struct file.
[filename, pathname] = uigetfile('*.mat', 'Choose the processed dvt or scored file. (*.mat)');
dvtFileName = fullfile(pathname, filename);
% !!!!!! NOTE !!!!!!! Put in error catch in case they hit cancel here.
% Load dvt file contents.
indRecStruct = load(dvtFileName).indRecStruct;

% Initialize variables - many are used in local and nested functions.    
if isfield(indRecStruct, 'events')
    events = indRecStruct.events;
    eventCount = length(events);
    pixelDvt = indRecStruct.pixelDVT;
    
    %Repainting pixelDvt.
    pathList = unique(events(:,3));
    pathList = pathList(abs(pathList)>0 & abs(pathList)<99);
    nPaths = length(pathList);
    for iPath = 1:nPaths
        runStarts = find(events(:,3) == pathList(iPath));
        for iRun = 1:length(runStarts)
            pixelDvt(events(runStarts(iRun),1):events(runStarts(iRun)+1,1),end) = ...
                pathList(iPath);
        end
    end
else
    % Initialize
    eventCount = 0;
    events = zeros(0,3);
    pixelDvt = indRecStruct.processedDVT;
    pixelDvt(:,end+1) = 99; %The last row will be updated with run label values.
end

% Save original pointer type to be used in mouse hovering
originalPointerType = get(gcf,'Pointer');
% Set function to run on close.
set(overallWindow, 'DeleteFcn', @hCloseFunction);
% Set event list print format.
format short g;
% Initialize variables so they exist across functions.
runStartClick = []; %Saves the click X location for graphing.
runEndClick = []; %Saves the click X location for graphing.
eventBuffer = [];
whichPath = 99; % Tag value for the path to plot when only plotting one.
% Sets the mouse back to "ignore-clicks-on-the-grid" mode.
stateResetter;
buttonResetter;

%STUB LOAD MOVIE HERE
% [filename, pathname] = uigetfile('*.avi', 'Choose the video file.');
% handles.aviFileName = fullfile(pathname, filename);
% !!!!!! NOTE !!!!!!! Put in error catch in case the hit cancel here.

% Change units to normalized so components resize automatically.
set([overallWindow, hSeeEventList, hSaveEventList, hMarkRun, hCleanRun,...
    hDirtyRun, hIndecisionPoint, hSelectAndRelabelRun,...
    hEraseLast, hSelectAndEraseRun, hSelectAndEraseEvent...
    hLightChooser, hLight1, hLight2, hLight3,hLight4,...
    hTimeSlider, hWindowTimeSizeLabel, hWindowTimeSize, hAddedRunLabels,...
    hAxes, hPlotOnePathOnly, hRunChooser, hAbort,hNoReward...
    hRun1,hRun2,hRun3,hRun4,hRun5,hRun6,hRun7,hRun8,hRun9,hRun10],...
    'Units','normalized');

% Assign the GUI a name to appear in the window title.
set(overallWindow,'Name','TT Position Tracker');
% Move the GUI to the center of the screen.
movegui(overallWindow,'center');
%Draw first set of data using default values.
drawPoints;
% Make the GUI visible.
guiReady = true;
set(overallWindow,'Visible','on');


%%  Callbacks for posTracker2014
% Callbacks automatically have access to variables nested at a lower level.
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
        
        stateResetter;
        buttonResetter;
        
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

        disableRunButtons;
        set(hNoReward,'Enable','off');
        set(hPlotOnePathOnly,'Enable','on');

        drawPoints;
    end

    function seeEventList_Callback(source,eventdata)
        % Print out spiral_events list to the command prompt.
        display(events);
        disp([num2str(eventCount),' lines right now.']);
    end

    function saveEventList_Callback(source,eventdata)
        % Save spiral_events labeled and the pixelDvt matrix (painted dvt).
        indRecStruct.pixelDVT = pixelDvt;
        indRecStruct.events = events;
        uisave('indRecStruct',strcat(filename(1:30),'_Events'));
    end

    function markRun_Callback(source,eventdata)
        % Enable run labeling. Click twice on grid after this.
        clickStateEnable;
        drawPoints;
    end

    function cleanRun_Callback(source,eventdata)
        % Mark current run as clean.
        set(hCleanRun,'Tag','1');
        set(hDirtyRun,'Value',0);
        
        % if in plot one path only mode
        if get(hPlotOnePathOnly,'Value')
            set(hSelectAndRelabelRun,'Enable','on');
            set(hEraseLast,'Enable','on');
            set(hSelectAndEraseRun,'Enable','on');
            set(hSelectAndEraseEvent,'Enable','on');
        else
            enableRunButtons;
        end      
        
        % if is deselected, undo all of the above
        if ~get(hCleanRun,'Value')
            if get(hPlotOnePathOnly,'Value')
                set(hSelectAndRelabelRun,'Enable','off');
                set(hEraseLast,'Enable','off');
                set(hSelectAndEraseRun,'Enable','off');
                set(hSelectAndEraseEvent,'Enable','off');
            else
                disableRunButtons;
            end
        end
        drawPoints;
    end

    function dirtyRun_Callback(source,eventdata)
        % Mark current run as dirty.
        set(hCleanRun,'Tag','-1');
        set(hCleanRun,'Value',0);

        % if in plot one path only mode
        if get(hPlotOnePathOnly,'Value')
            set(hSelectAndRelabelRun,'Enable','on');
            set(hEraseLast,'Enable','on');
            set(hSelectAndEraseRun,'Enable','on');
            set(hSelectAndEraseEvent,'Enable','on');
        else
            enableRunButtons;
        end
                
        % if is deselected, undo all of the above
        if ~get(hDirtyRun,'Value')
            if get(hPlotOnePathOnly,'Value')
                set(hSelectAndRelabelRun,'Enable','off');
                set(hEraseLast,'Enable','off');
                set(hSelectAndEraseRun,'Enable','off');
                set(hSelectAndEraseEvent,'Enable','off');
            else
                disableRunButtons;
            end
        end
        drawPoints;
    end

    function indecisionPoint_Callback(source,eventdata)
        % Enable indecision point labeling. Click once on grid after this.
        clickStateEnable;
        drawPoints;
    end

    function selectAndRelabelRun_Callback(source,eventdata)
         set(gcf,'WindowButtonMotionFcn',@detectMouseHover);
         
         set(hMarkRun,'Enable','off');
         set(hIndecisionPoint,'Enable','off');
         set(hSelectAndRelabelRun,'Enable','off');
         set(hEraseLast,'Enable','off');
         set(hSelectAndEraseRun,'Enable','off');
         set(hSelectAndEraseEvent,'Enable','off');
         
         drawPoints;
     end

    function selectAndEraseRun_Callback(source,eventdata)
         set(gcf,'WindowButtonMotionFcn',@detectMouseHover);
         
         set(hMarkRun,'Enable','off');
         set(hIndecisionPoint,'Enable','off');
         set(hSelectAndRelabelRun,'Enable','off');
         set(hEraseLast,'Enable','off');
         set(hSelectAndEraseRun,'Enable','off');
         set(hSelectAndEraseEvent,'Enable','off');
         
         drawPoints;
    end

    function selectAndEraseEvent_Callback(source,eventdata)
        set(gcf,'WindowButtonMotionFcn',@detectMouseHover);
        
        set(hMarkRun,'Enable','off');
        set(hIndecisionPoint,'Enable','off');
        set(hSelectAndRelabelRun,'Enable','off');
        set(hEraseLast,'Enable','off');
        set(hSelectAndEraseRun,'Enable','off');
        set(hSelectAndEraseEvent,'Enable','off');
        
        drawPoints;
    end

    function eraseLastEvent_Callback(source,eventdata)
        if eventCount <=1
            eventCount = 0;
            events = zeros(0,3);
            pixelDvt(:,end) = 99;
        else
            lastLabel = events(end,end);
            if lastLabel == 777
                %Trash 3 lines.
                eventCount = eventCount-3;
                pixelDvt(events(end-2,1):events(end-1,1),end) = 99;
                events = events(1:end-3,:);
            elseif abs(lastLabel) > 100
                %Trash 2 lines
                eventCount = eventCount-2;
                pixelDvt(events(end-1,1):events(end,1),end) = 99;
                events = events(1:end-2,:);
            else
                eventCount = eventCount-1;
                events = events(1:end-1,:);
            end
        end
        
        % Reset and get ready for another event.
        stateResetter;
        buttonResetter;
        drawPoints;
    end

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

    function addedRunLabels_Callback(source,eventdata)
        % Toggle the run buttons to tag runs either 1-10 or 11-20.
        if get(source,'value')
            set(hRun1,'Tag','11');
            set(hRun1,'String','Run 11');
            
            set(hRun2,'Tag','12');
            set(hRun2,'String','Run 12');
            
            set(hRun3,'Tag','13');
            set(hRun3,'String','Run 13');
            
            set(hRun4,'Tag','14');
            set(hRun4,'String','Run 14');
            
            set(hRun5,'Tag','15');
            set(hRun5,'String','Run 15');
            
            set(hRun6,'Tag','16');
            set(hRun6,'String','Run 16');
            
            set(hRun7,'Tag','17');
            set(hRun7,'String','Run 17');
            
            set(hRun8,'Tag','18');
            set(hRun8,'String','Run 18');
            
            set(hRun9,'Tag','19');
            set(hRun9,'String','Run 19');
            
            set(hRun10,'Tag','20');
            set(hRun10,'String','Run 20');
        else
            set(hRun1,'Tag','1');
            set(hRun1,'String','Run 1');
            
            set(hRun2,'Tag','2');
            set(hRun2,'String','Run 2');
            
            set(hRun3,'Tag','3');
            set(hRun3,'String','Run 3');
            
            set(hRun4,'Tag','4');
            set(hRun4,'String','Run 4');
            
            set(hRun5,'Tag','5');
            set(hRun5,'String','Run 5');
            
            set(hRun6,'Tag','6');
            set(hRun6,'String','Run 6');
            
            set(hRun7,'Tag','7');
            set(hRun7,'String','Run 7');
            
            set(hRun8,'Tag','8');
            set(hRun8,'String','Run 8');
            
            set(hRun9,'Tag','9');
            set(hRun9,'String','Run 9');
            
            set(hRun10,'Tag','10');
            set(hRun10,'String','Run 10');
        end
    end

    function plotOnePathOnly_Callback(source,eventdata)
        % Plot only one type of runs.
        whichPath = 99;
        set(hMarkRun,'Enable','off');
        set(hCleanRun,'Enable','on');
        set(hDirtyRun,'Enable','on');
        set(hIndecisionPoint,'Enable','off');
        set(hSelectAndRelabelRun,'Enable','off');
        set(hEraseLast,'Enable','off');
        set(hSelectAndEraseRun,'Enable','off');
        set(hSelectAndEraseEvent,'Enable','off');
        
        enableRunButtons;
        drawPoints;
        
        % Reset run buttons if unselected.
        if ~get(hPlotOnePathOnly,'Value')
            set(hRunChooser,'SelectedObject',[]);
            buttonResetter;
        end
        drawPoints;
    end

    function runChooser_Callback(source,eventdata) %#ok<*INUSL>
        % Runs when a run selection button is pressed. Does a lot of work
        % with recording events.
        
        if get(hSelectAndRelabelRun,'Value')
            % Find which button is pressed.
            hRunChoice = get(hRunChooser,'SelectedObject');
            runChoice = str2double(get(hRunChoice,'Tag'))*str2double(get(hCleanRun,'Tag'));
            % Calculate the run markers.
            nEventsToAdd = size(eventBuffer,1);             %#ok<*NODEF>
            eventBuffer(1:2,3) = [runChoice;runChoice+sign(runChoice)*100];            
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
            % Reset and get ready for another event.
            whichPathTmp = whichPath;
            stateResetter;
            buttonResetter;
            whichPath = whichPathTmp;
            set(hSelectAndRelabelRun,'Value',0);
        elseif get(hPlotOnePathOnly,'Value')
            % If for the plot one path selection.
            set(hSelectAndRelabelRun,'Enable','on');
            set(hSelectAndEraseRun,'Enable','on');
            set(hEraseLast,'Enable','on');
            whichPath = str2double(get(eventdata.NewValue,'Tag'));
%         elseif eventdata.NewValue == hNoReward
%             % Reset and get ready for another event.
%             stateResetter;
%             buttonResetter;
        else
            %Inactivate appropriate buttons.
            set(hCleanRun,'Enable','off');
            set(hDirtyRun,'Enable','off');
            
            % Find which button is pressed.
            hRunChoice = get(hRunChooser,'Selectedobject');
            runChoice = str2double(get(hRunChoice,'Tag'))*str2double(get(hCleanRun,'Tag'));
            
            % Calculate the run markers.
            eventBuffer(:,3) = [runChoice;runChoice+sign(runChoice)*100];
            
            % Check to make sure this exact run has not already been added. I
            % do this this way because now if they push a different run button
            % it just overwrites the previous button - a simple way for users
            % to deal with a misclick.
            if size(events,1) >= 2
                if ~(eventBuffer (:,1:2) == events(end-1:end,1:2))
                    eventCount = eventCount + 2;
                end
            else
                eventCount = eventCount + 2;
            end
            % Add the run markers to the spiral_events list.
            events(eventCount-1:eventCount,:) = eventBuffer;
            
            % Paint the tracking file appropriately.
            pixelDvt(events(end-1,1):events(end,1),end) = ...
                events(end-1,3);
            
%             % REWARD CODE - Commented out.
%             % Get one more click for the reward location or click the "No Reward" button.
%             % set(hNoReward,'Enable','on');
% 
%             % Start Mouse Hovering Detection
%             % Asks the user to position the cursor with the mouse, and then returns
%             % those coordinates upon any click or key press.
%             set(gcf,'WindowButtonMotionFcn',@detectMouseHover);
            stateResetter;
            buttonResetter;
        end
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
        % Handles mouse clicks - this it the function that does all of the
        % work - any time a click is recorded the state of the buttons must
        % be determined and the appropriate action taken.
        
        % Get the location of the mouse at time of the click.
        clickLocation = get(hAxes,'CurrentPoint');
        
        % Get info about the status of the data in the interface.
        [firstPoint, lastPoint, xLightCol, yLightCol] = statusSet();
        
        % Typical event labeling step to determine which click just
        % happened. Doesn't happen in selectAndEraseRun.
        if ~get(hSelectAndEraseRun,'Value') && ~get(hSelectAndRelabelRun,'Value')
            nClicks = str2double(get(hAxes,'Tag'));
            nClicks = nClicks+1;
            set(hAxes,'Tag',num2str(nClicks));
        end
        
        % Determine plotted points.
        allLinesPlotted = [];
        if get(hSelectAndEraseRun,'Value') || get(hSelectAndRelabelRun,'Value')
            % Plotted points are only certain runs.
            
            relevantTypes = false(size(events,1),1); %#ok<*PREALL> % Initialize
            % Gives all relevant spiral_events for the chosen time window.
            relevantTimes = (firstPoint <= events(:,1)) &...
                (events(:,1) <= lastPoint);
            
            if get(hCleanRun,'Value') == 1
                % Gives all clean run spiral event starts.
                relevantTypes = (0 < events(:,3)) &...
                    (events(:,3) < 100);
            elseif get(hDirtyRun,'Value') == 1
                % Gives all dirty run spiral event starts.
                relevantTypes = (0 > events(:,3)) &...
                    (events(:,3) > -100);
            else
                relevantTypes = (0 < abs(events(:,3))) &...
                    (abs(events(:,3)) < 100);
            end
            if ~isempty(get(hRunChooser,'SelectedObject'))
                % Gives all runs of selected run type.
                relevantTypes = relevantTypes &...
                    abs(events(:,3)) == whichPath;
            end
            % Get the starts and ends.
            relevantStartEvents = find(relevantTimes & relevantTypes);
            relevantEndEvents = find([0;relevantTypes(1:end-1)] & relevantTimes);
            
            % No runs being plotted currently
            if isempty(relevantStartEvents) && isempty(relevantEndEvents)
                uiwait(msgbox(...
                    {'There are no runs being plotted to be selected!',...
                    'Please make sure the run you want to effect is plotted.'},...
                    'What run are you selecting?!','warn'));
                abort_Callback;
                return;
            end
            
            % Include a partial run at the start of the time window being plotted.
            if min(relevantEndEvents) < min(relevantStartEvents)
                allLinesPlotted = [allLinesPlotted;...
                    (firstPoint:events(relevantEndEvents(1),1))'];
            end
            
            for iEvent = 1:length(relevantEndEvents);
                allLinesPlotted = [allLinesPlotted;...
                    (events(relevantEndEvents(iEvent)-1,1):...
                    events(relevantEndEvents(iEvent),1))']; %#ok<*AGROW>
            end
            
            % Include a partial run at the end of the time window being plotted.
            if isempty(iEvent) 
                iEvent = 1;
                allLinesPlotted = [allLinesPlotted;...
                    (events(relevantStartEvents(iEvent),1):lastPoint)'];
            elseif iEvent < length(relevantStartEvents)
                iEvent = iEvent+1;
                allLinesPlotted = [allLinesPlotted;...
                    (events(relevantStartEvents(iEvent),1):lastPoint)'];
            end
            
        elseif get(hSelectAndEraseEvent,'Value')
            % For the indecision points, at least for now.
            
            relevantTypes = events(:,3) == 0; %#ok<*PREALL> % Initialize
            % Gives all relevant spiral_events for the chosen time window.
            relevantTimes = (firstPoint <= events(:,1)) &...
                (events(:,1) <= lastPoint);
            relevantEvents = find(relevantTimes & relevantTypes);
            allLinesPlotted = [events(relevantEvents,1)]';
            if isempty(relevantEvents)
                uiwait(msgbox(...
                    {'There are no events being plotted to be selected!',...
                    'Please make sure the event you want to effect is plotted.'},...
                    'What event are you selecting?!','warn'));
%                 uiwait(msgbox(...
%                     {'Are you being serious right now??!'},...
%                     'Get it together.','warn'));
                abort_Callback;
                return;
            end
        else
            % Normal Situation
            allLinesPlotted = [firstPoint:lastPoint]'; %#ok<*NBRAK>           
        end
        
        % Determine point currently plotted that click is closest to.
        runDist = dist(clickLocation(1,1:2),pixelDvt(allLinesPlotted,[xLightCol,yLightCol])');
        [~, ind] = min(runDist);
        % this is the row index of the user-selected point, with
        % respect to the full 'pos' matrix.
        dvtRow = allLinesPlotted(ind(1));
                
        
        % Process click. Two general cases - selecting a run or marking an
        % event.
        if get(hSelectAndEraseRun,'Value') || get(hSelectAndRelabelRun,'Value')
            % The selectAndEraseRun & selectAndRelabelRun code
            
            % selectAndRelabelRun code - identify the lines, blast them,
            % store a copy in eventBuffer. Then let the user select
            % good/bad and run #. Only change is not deleting event buffer.
            % Store bad data there.
            
            % Identify the appropriate run.
            [startDiff,indStartEvents] = min(abs(dvtRow-events(relevantStartEvents,1)));
            [endDiff,indEndEvents] = min(abs(dvtRow-events(relevantEndEvents,1)));
            
            % Finding the relevant rows in the position file & spiral_events.
            if isempty(endDiff) | (startDiff <= endDiff)
                % Click was closer to the start of the run. Grab the run.
                badEventRowInSpiral = relevantStartEvents(indStartEvents);
            else
                % Click was closer to the end of the run. Grab the run.
                badEventRowInSpiral = relevantEndEvents(indEndEvents)-1;
            end
            blastRunStartDvtRow = events(badEventRowInSpiral,1);
            blastRunEndDvtRow = events(badEventRowInSpiral+1,1);
            
            % Repainting the points in pixelDvt to unlabeled.
            pixelDvt(blastRunStartDvtRow:blastRunEndDvtRow,end) = 99;

            % Erase the old event. Will save the run in eventBuffer and add 
            % it to the end of spiral_events if we are relabeling.
            if eventCount == 2
                eventCount = 0;
                eventBuffer = events;
                drawPoints;
                events = zeros(0,3);
            elseif size(events,1) >= badEventRowInSpiral+2
                if events(badEventRowInSpiral+2,end) == 777; % Remove reward label too.
                eventCount = eventCount-3;
                eventBuffer = events(badEventRowInSpiral:badEventRowInSpiral+2,:);
                drawPoints;
                events = events([1:badEventRowInSpiral-1,badEventRowInSpiral+3:end],:);
                else
                    eventCount = eventCount-2;
                    eventBuffer = events(badEventRowInSpiral:badEventRowInSpiral+1,:);
                    drawPoints;
                    events = events([1:badEventRowInSpiral-1,badEventRowInSpiral+2:end],:);
                end
            else
                eventCount = eventCount-2;
                eventBuffer = events(badEventRowInSpiral:badEventRowInSpiral+1,:);
                events = events([1:badEventRowInSpiral-1,badEventRowInSpiral+2:end],:);
            end
       
            % If we are running select and relabel, enable correct buttons.
            if get(hSelectAndRelabelRun,'Value')
                set(hCleanRun,'Enable','on');
                set(hDirtyRun,'Enable','on');
                
                % Reset mouse click portion but not complete state reset.
                set(hAxes,'Tag','-1');
                set(overallWindow,'WindowButtonMotionFcn',''); % Hope this undoes the callback.
                set(gcf,'Pointer',originalPointerType);
                % If doing nothing, you need to tell it to do nothing, or else
                % the other condition for 'WindowButtonDownFcn' will be exceuted
                %  from the previous mouse position and that is not desirable.
                set(gcf,'WindowButtonDownFcn','');
                set(hRunChooser,'SelectedObject',[]);
                set(hPlotOnePathOnly,'Value',0);
                
                disableRunButtons;
                set(hNoReward,'Enable','off');
            else
                stateResetter; % Clear the event out, we don't need it.
                buttonResetter;
                drawPoints;
            end
        
        elseif get(hSelectAndEraseEvent,'Value')
            % Indecision point label removal
            
            % Identify the appropriate event.
            [~,indEvents] = min(abs(dvtRow-events(relevantEvents,1)));
            
            % Finding the relevant rows in the position file & spiral_events.
            badEventRowInSpiral = relevantEvents(indEvents);
            
            % Erase the old event.
            if eventCount == 1
                eventCount = 0;
                events = zeros(0,3);
            else
                eventCount = eventCount-1;
                events = events([1:badEventRowInSpiral-1,badEventRowInSpiral+1:end],:);
            end
            
            stateResetter; % Clear the event out, we don't need it.
            buttonResetter;
            drawPoints;
            
        else
            % Typical event labeling step to determine which click just
            % happened.
            if nClicks == 1
                % Indecision point or run?
                if get(hIndecisionPoint,'Value')
                    % Record indecision location, reset and get ready for another event.
                    eventCount = eventCount+1;
                    events(eventCount,:) = [dvtRow, pixelDvt(dvtRow,2), 0];
                    
                    % Sets the mouse back to "ignore-clicks-on-the-grid" mode.
                    stateResetter;
                    buttonResetter;
                    
                else %Normal scenario of labeling a run
                    eventBuffer = [dvtRow, pixelDvt(dvtRow,2)];
                    runStartClick = clickLocation;
                end
            elseif nClicks == 2
                if dvtRow < eventBuffer(1)
                    % Clicked backwards.
                    % Pop up a message box to explain and reset everything.
                    uiwait(msgbox(...
                        {'Your second click happend before your first click in time!',...
                        'Maybe he is running the other direction?',...
                        'This event will be cleared. Please try again.'},...
                        'How is he running backwards in time?!','warn'));
                    % Sets the mouse back to "ignore-clicks-on-the-grid" mode.
                    stateResetter;
                    buttonResetter;
                else
                    % Normal situation
                    eventBuffer(2,:) = [dvtRow, pixelDvt(dvtRow,2)];
                    runEndClick = clickLocation;
                    
                    % Sets the mouse back to "ignore-clicks-on-the-grid" mode.
                    set(overallWindow,'WindowButtonMotionFcn',''); % Hope this undoes the callback.
                    set(gcf,'Pointer',originalPointerType);
                    % If doing nothing, you need to tell it to do nothing, or else
                    % the other condition for 'WindowButtonDownFcn' will be exceuted
                    %  from the previous mouse position and that is not desirable.
                    set(gcf,'WindowButtonDownFcn','');
                    set(hRunChooser,'SelectedObject',[]);
                    
                    % Activate appropriate buttons.
                    set(hCleanRun,'Enable','on');
                    set(hDirtyRun,'Enable','on');
                end
            elseif nClicks == 3
                % Record reward location, reset and get ready for another event.
                eventCount = eventCount+1;
                events(eventCount,:) = [dvtRow, pixelDvt(dvtRow,2), 777];
                
                % Sets the mouse back to "ignore-clicks-on-the-grid" mode.
                stateResetter;
                buttonResetter;
            end
            drawPoints;
        end

    end

%% GUI close callback
    function hCloseFunction(handles,object,eventdata)
        % Save spiral_events labeled and the pixelDvt matrix (painted dvt).
%         save('posTrackerTMPSAVEFILE', 'events', 'pixelDvt');
    end

%%  Utility functions for posTracker2014
    function newAxesPos = estimateNewPosition()
        %% This function is used to get the new position of the uipanel, to be used for GUIs that are Resizable.
        
        guiPosition = get(overallWindow,'Position');
        axesPosition = get(hAxes,'Position');
        
        newAxesPos(1) = guiPosition(3)*axesPosition(1);
        newAxesPos(2) = guiPosition(4)*axesPosition(2);
        newAxesPos(3) = guiPosition(3)*axesPosition(3);
        newAxesPos(4) = guiPosition(4)*axesPosition(4);
        
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

    function clickStateEnable()
        % Start mouse hovering detection & turn off appropriate buttons.
        % asks the user to position the cursor with the mouse, and then returns
        % those coordinates upon any click or key press.
        set(hAxes,'Tag','0');
        set(gcf,'WindowButtonMotionFcn',@detectMouseHover);
        set(hMarkRun,'Enable','off');
        set(hIndecisionPoint,'Enable','off');
        set(hSelectAndRelabelRun,'Enable','off');
        set(hEraseLast,'Enable','off');
        set(hSelectAndEraseRun,'Enable','off');
        set(hPlotOnePathOnly,'Enable','off');
    end

    function stateResetter()
        % Reset and get ready for another event.
        % Sets the mouse back to "ignore-clicks-on-the-grid" mode.
        set(hAxes,'Tag','-1');
        set(overallWindow,'WindowButtonMotionFcn',''); % Hope this undoes the callback.
        set(gcf,'Pointer',originalPointerType);
        % If doing nothing, you need to tell it to do nothing, or else
        % the other condition for 'WindowButtonDownFcn' will be exceuted
        %  from the previous mouse position and that is not desirable.
        set(gcf,'WindowButtonDownFcn','');
        set(hRunChooser,'SelectedObject',[]);
        set(hPlotOnePathOnly,'Value',0);
        runStartClick = []; %Saves the click X location for graphing.
        runEndClick = []; %Saves the click X location for graphing.
        eventBuffer = [];
        whichPath = 99; % Tag value for the path to plot when only plotting one.
    end

    function buttonResetter()
        set(hCleanRun,'Tag','0');
        
        set(hMarkRun,'Value',0);
        set(hCleanRun,'Value',0);
        set(hDirtyRun,'Value',0);
        set(hIndecisionPoint,'Value',0);     
        
        set(hMarkRun,'Enable','on');
        set(hCleanRun,'Enable','off');
        set(hDirtyRun,'Enable','off');
        set(hIndecisionPoint,'Enable','on');
        
        disableRunButtons;
        set(hNoReward,'Enable','off');

        set(hSelectAndRelabelRun,'Value',0);
        set(hSelectAndRelabelRun,'Enable','on');
        
        set(hSelectAndEraseRun,'Value',0);
        set(hSelectAndEraseRun,'Enable','on');
        set(hSelectAndEraseEvent,'Value',0);
        set(hSelectAndEraseEvent,'Enable','on');
        
        set(hEraseLast,'Enable','on');
        
        set(hPlotOnePathOnly,'Enable','on');
    end

    function enableRunButtons()
        %Enable run buttons.
        set(hRun1,'Enable','on');
        set(hRun2,'Enable','on');
        set(hRun3,'Enable','on');
        set(hRun4,'Enable','on');
        set(hRun5,'Enable','on');
        set(hRun6,'Enable','on');
        set(hRun7,'Enable','on');
        set(hRun8,'Enable','on');
        set(hRun9,'Enable','on');
        set(hRun10,'Enable','on');
    end

    function disableRunButtons()
        %Enable run buttons.
        set(hRun1,'Enable','off');
        set(hRun2,'Enable','off');
        set(hRun3,'Enable','off');
        set(hRun4,'Enable','off');
        set(hRun5,'Enable','off');
        set(hRun6,'Enable','off');
        set(hRun7,'Enable','off');
        set(hRun8,'Enable','off');
        set(hRun9,'Enable','off');
        set(hRun10,'Enable','off');
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
        
        % For plotting only one path type
        if get(hPlotOnePathOnly,'Value')
            if isempty(get(hRunChooser,'SelectedObject'))
                % Plot all run types.
                if get(hCleanRun,'Value') == 1
                    %All clean runs of any type.
                    onePath = find(pixelDvt(firstPoint+changeStarts,end)>0 &...
                        pixelDvt(firstPoint+changeStarts,end)<99); % Clean
                elseif get(hDirtyRun,'Value') == 1
                    %All dirty runs of any type.
                    onePath = find(pixelDvt(firstPoint+changeStarts,end)<0 &...
                        pixelDvt(firstPoint+changeStarts,end)>-99); % Dirty
                else
                    % All unlabeled runs.
                    onePath = find(pixelDvt(firstPoint+changeStarts,end) == whichPath); %Unlabeled
                end
            else
                if get(hCleanRun,'Value') == 1
                    onePath = find(pixelDvt(firstPoint+changeStarts,end) == whichPath); % Clean one run
                elseif get(hDirtyRun,'Value') == 1
                    onePath = find(pixelDvt(firstPoint+changeStarts,end) == -1*whichPath); % Dirty one run
                else
                    onePath = find(abs(pixelDvt(firstPoint+changeStarts,end)) == whichPath); %Unlabeled
                end
            end
            changeStarts = changeStarts(onePath);
            changeEnds = changeEnds(onePath);
            
            % Prompt to use run buttons.
            if ~(get(hSelectAndEraseEvent,'Value') ||...
                    get(hSelectAndEraseRun,'Value') ||...
                    get(hSelectAndRelabelRun,'Value'))
                text(0,yAxisRange+10,'Click the button for the path you want to see.',...
                    'FontSize',15,'FontWeight','bold','Color','k');
            end
        end
        
        % Plots the labeled portions in the correct colors.
        for iSegment = 1:length(changeStarts)
            plot(pixelDvt(firstPoint+changeStarts(iSegment):firstPoint+changeEnds(iSegment),xLightCol),...
                pixelDvt(firstPoint+changeStarts(iSegment):firstPoint+changeEnds(iSegment),yLightCol),'.',...
                'Color',pen(abs(pixelDvt(firstPoint+changeStarts(iSegment),end)),:));
        end
        
        % Plot reward locations.
        events = indRecStruct.events;
        rewardLocations = find(events(:,end) == 777); % Reward spots
        onScreenRewardLocs = (firstPoint<= events(rewardLocations,1)) &...
            (events(rewardLocations,1) <= lastPoint);
        rewardLocsToPlot = rewardLocations(onScreenRewardLocs);
        for iReward = 1:length(rewardLocsToPlot)
            plot(pixelDvt(events(rewardLocsToPlot(iReward),1),xLightCol),...
                pixelDvt(events(rewardLocsToPlot(iReward),1),yLightCol),...
                'v','Color',[0.9,0,0.9]);
        end
        
        % Plot indecision locations.
        indecisionLocations = find(events(:,end) == 0); % Reward spots
        onScreenIndecisionLocs = (firstPoint<= events(indecisionLocations,1)) &...
            (events(indecisionLocations,1) <= lastPoint);
        indecisionLocsToPlot = indecisionLocations(onScreenIndecisionLocs);
        for iIndecision = 1:length(indecisionLocsToPlot)
            plot(pixelDvt(events(indecisionLocsToPlot(iIndecision),1),xLightCol),...
                pixelDvt(events(indecisionLocsToPlot(iIndecision),1),yLightCol),...
                's','Color',[0.9,0,0.9]);
        end
        
        % Plots the start and end points for the time segment.
        plot(pixelDvt(firstPoint,xLightCol),...
            pixelDvt(firstPoint,yLightCol),'g+','linewidth',3);
        plot(pixelDvt(lastPoint,xLightCol),...
            pixelDvt(lastPoint,yLightCol),'r+','linewidth',3);
        set(hAxes,'Tag',num2str(nClicks)); %Deals with weird fact that the plot fn resets axes 'Tag' value to ''

        % Plot conditional instruction text.
        if get(hSelectAndRelabelRun,'Value') && abs(str2double(get(hCleanRun,'Tag')))
            text(0,yAxisRange+10,'Click the path number.',...
                'FontSize',15,'FontWeight','bold','Color','k');
        elseif get(hSelectAndRelabelRun,'Value') && ~isempty(eventBuffer)
            text(0,yAxisRange+10,'Choose if it is a clean or dirty run to the right.',...
                'FontSize',15,'FontWeight','bold','Color','k');
        elseif get(hSelectAndRelabelRun,'Value')
            text(0,yAxisRange+10,'Click the run you want to relabel.',...
                'FontSize',15,'FontWeight','bold','Color','k');
        elseif get(hSelectAndEraseRun,'Value')
                text(0,yAxisRange+10,'Click the run you want to erase.',...
                    'FontSize',15,'FontWeight','bold','Color','k');
        elseif get(hSelectAndEraseEvent,'Value')
                    text(0,yAxisRange+10,'Click the event you want to erase.',...
                        'FontSize',15,'FontWeight','bold','Color','k');
        end
        
        % Plot conditionally based on the current number of clicks.
        switch nClicks
            case 0
                if get(hIndecisionPoint,'Value')
                    text(0,yAxisRange+10,'Click the indecision point.',...
                        'FontSize',15,'FontWeight','bold','Color','k');                  
                else
                    text(0,yAxisRange+10,'Click the run start.',...
                        'FontSize',15,'FontWeight','bold','Color','k');
                end
            case 1
                text(0,yAxisRange+10,'Click the run end',...
                    'FontSize',15,'FontWeight','bold','Color','k');
                % Print a plus on the clicked point, to be used as an indicator
                text(runStartClick(1,1),runStartClick(1,2),'X','FontSize',15,'FontWeight','bold','Color','g');
            case 2
                if ~isempty(get(hRunChooser,'Selectedobject'))
                    text(0,yAxisRange+10,'Revise your run choice or select reward location.',...
                        'FontSize',15,'FontWeight','bold','Color','k');
                    % Print pluses on the clicked points, to be used as an indicator
                    text(runStartClick(1,1),runStartClick(1,2),'X','FontSize',15,'FontWeight','bold','Color','g');
                    text(runEndClick(1,1),runEndClick(1,2),'X','FontSize',15,'FontWeight','bold','Color','r');
                elseif str2double(get(hCleanRun,'Tag')) == 0
                    text(0,yAxisRange+10,'Choose if it is a clean or dirty run to the right.',...
                        'FontSize',15,'FontWeight','bold','Color','k');
                    % Print pluses on the clicked points, to be used as an indicator
                    text(runStartClick(1,1),runStartClick(1,2),'X','FontSize',15,'FontWeight','bold','Color','g');
                    text(runEndClick(1,1),runEndClick(1,2),'X','FontSize',15,'FontWeight','bold','Color','r');
                else
                    text(0,yAxisRange+10,'Choose the path number.',...
                        'FontSize',15,'FontWeight','bold','Color','k');
                    % Print pluses on the clicked points, to be used as an indicator
                    text(runStartClick(1,1),runStartClick(1,2),'X','FontSize',15,'FontWeight','bold','Color','g');
                    text(runEndClick(1,1),runEndClick(1,2),'X','FontSize',15,'FontWeight','bold','Color','r');
                end
        end
    end
end