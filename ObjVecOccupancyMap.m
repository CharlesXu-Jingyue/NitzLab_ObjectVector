%% Object Vector Occupancy Map
% Charles Xu @ UCSD, v1, 20221114
% Create Occupancy Map in allocentric coordinates (room or
% object-centered)

%% Load processed indRecStruct file
clear

[matFileName, matPathName] = uigetfile('*.mat', 'Choose the indRecStruct file.');
load(fullfile(matPathName,matFileName))

dvtRoom = indRecStruct.processedDVT(:,9:10);
dvtRoom(dvtRoom==1) = NaN;
dvtObj = indRecStruct.objVec.processedDVT(:,9:10);
dvtObj(dvtObj==1) = NaN;

width = 600; % specify resolution
height = 600;
gridSize = 100; % specify grid size

%% Create occupancy map (Update to get occupancy map for each inner run)
% Occupancy matrix in room coordinates
hRoom = hist3(dvtRoom,'Edges',{(1:gridSize:width)-1,(1:gridSize:height)-1});
% Plot occupancy map
hist3(dvtRoom,'Edges',{(1:gridSize:width)-1,(1:gridSize:height)-1},'CDataMode','auto')
ylabel('Door Side')
colorbar
view(2)

% Occupancy matrix in object-centered coordinates
hObj = hist3(dvtObj,'Edges',{(1:gridSize:width)-1,(1:gridSize:height)-1});
% Plot occupancy map
hist3(dvtObj,'Edges',{(1:gridSize:width)-1,(1:gridSize:height)-1},'CDataMode','auto')
ylabel('Door Side')
colorbar
view(2)

% mapRoomCentered = occupancyMap(width,height);
% mapObjCentered = occupancyMap(2*width,2*height);
% 
% [X,Y] = meshgrid(1:width, 1:height);
% setOccupancy(mapRoomCentered, [X(:) Y(:)], hRoom(:))
% setOccupancy(mapObjCentered, [X(:) Y(:)], hObj(:))

