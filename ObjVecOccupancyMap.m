%% Object Vector Occupancy Map
% Charles Xu @ UCSD, v1, 20221114
% Create Occupancy Map in allocentric coordinates (room or
% object-centered)

%% Load processed indRecStruct file
clear

[matFileName, matPathName] = uigetfile('*.mat', 'Choose the indRecStruct file.');
load(fullfile(matPathName,matFileName))

%% Initialize
dvtRoom = indRecStruct.processedDVT(:,[1,2,9,10]);
dvtRoom(dvtRoom==1) = NaN;
dvtRoom(1,1) = 1;
dvtObj = indRecStruct.objVec.processedDVT(:,[1,2,9,10]);
dvtObj(dvtObj==1) = NaN;
dvtObj(1,1) = 1;
objPosRoom = indRecStruct.objPosition;
objPosObj = indRecStruct.objVec.objPosition;
inner = indRecStruct.event.inner;

width = 600; % specify resolution
height = 600;
gridSize = 10; % specify grid size
runNum = 19; % specify event number to plot

%% Create occupancy map (Update to get occupancy map for each inner run)
if runNum <= size(inner,1)
    dvtRoomRun = dvtRoom(dvtRoom(:,2)>=inner(runNum,1) & dvtRoom(:,2)<=inner(runNum,2),:);
    % Occupancy matrix in room coordinates
    hRoom = hist3(dvtRoomRun(:,3:4),'Edges',{(1:gridSize:width)-1,(1:gridSize:height)-1});
    % Plot occupancy map
    figure
    hold on
    hist3(dvtRoomRun(:,3:4),'Edges',{(1:gridSize:width)-1,(1:gridSize:height)-1},'CDataMode','auto')
    ylabel('Door Side')
    colorbar
    view(2)
    hold off

    dvtObjRun = dvtObj(dvtObj(:,2)>=inner(runNum,1) & dvtObj(:,2)<=inner(runNum,2),:);
    % Occupancy matrix in object-centered coordinates
    hObj = hist3(dvtObjRun(:,3:4),'Edges',{(1:gridSize:width)-1,(1:gridSize:height)-1});
    % Plot occupancy map
    figure
    hold on
    hist3(dvtObjRun(:,3:4),'Edges',{(-width:gridSize:width)-1,(-height:gridSize:height)-1},'CDataMode','auto')
    ylabel('Door Side')
    colorbar
    view(2)
    hold off

    % mapRoomCentered = occupancyMap(width,height);
    % mapObjCentered = occupancyMap(2*width,2*height);
    % 
    % [X,Y] = meshgrid(1:width, 1:height);
    % setOccupancy(mapRoomCentered, [X(:) Y(:)], hRoom(:))
    % setOccupancy(mapObjCentered, [X(:) Y(:)], hObj(:))
end

