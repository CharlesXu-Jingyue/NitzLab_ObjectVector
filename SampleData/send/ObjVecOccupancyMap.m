%% ObjVecOccupancyMap
% Charles Xu @ UCSD, v1, 20221114
% ObjVecOccupancyMap creates an occupancy map in allocentric coordinates
% (world or object-centered) by reading from an event-processed indRecStruct

%% Load processed indRecStruct file
clear

[matFileName, matPathName] = uigetfile('*.mat', 'Choose the indRecStruct file.');
load(fullfile(matPathName,matFileName))

%% Initialize
dvtWorld = indRecStruct.world.processedDVT(:,[1,2,9,10]);
dvtWorld(dvtWorld==1) = NaN;
dvtWorld(1,1) = 1;
dvtObj = indRecStruct.objVec.processedDVT(:,[1,2,9,10]);
dvtObj(dvtObj==1) = NaN;
dvtObj(1,1) = 1;
objPosWorld = indRecStruct.world.objPosition;
objPosObj = indRecStruct.objVec.objPosition;
inner = indRecStruct.event.inner;

width = 600; % specify resolution
height = 600;
gridSize = 10; % specify grid size
runNum = 29; % specify event number to plot

%% Create occupancy map (Update to get occupancy map for each inner run)
if runNum <= size(inner,1)
    dvtWorldRun = dvtWorld(dvtWorld(:,2)>=inner(runNum,1) & dvtWorld(:,2)<=inner(runNum,2),:);
    dvtObjRun = dvtObj(dvtObj(:,2)>=inner(runNum,1) & dvtObj(:,2)<=inner(runNum,2),:);
    objPosWorldRun = objPosWorld(find(objPosWorld(:,2)>inner(runNum,1),1,'first'),:);
    objPosObjRun = objPosObj(find(objPosObj(:,2)>inner(runNum,1),1,'first'),:);
    
    % Raw trace in world-centered view
    figure
    hold on
    scatter(dvtWorldRun(:,3), dvtWorldRun(:,4), 'b')
    scatter(objPosWorldRun([3 5 7]), objPosWorldRun([4 6 8]), 'r', 'filled')
    xlim([0,width])
    ylim([0,height])
    hold off
    
    % Occupancy matrix in world coordinates
    hWorld = hist3(dvtWorldRun(:,3:4),'Edges',{(1:gridSize:width)-1,(1:gridSize:height)-1});
    % Plot occupancy map
    figure
    hold on
    hist3(dvtWorldRun(:,3:4),'Edges',{(1:gridSize:width)-1,(1:gridSize:height)-1},'CDataMode','auto')
    scatter(objPosWorldRun([3 5 7]), objPosWorldRun([4 6 8]), 'y', 'filled')
    ylabel('Door Side')
    colorbar
    title('World-centered')
    view(2)
    hold off

    % Occupancy matrix in object-centered coordinates
    hObj = hist3(dvtObjRun(:,3:4),'Edges',{(1:gridSize:width)-1,(1:gridSize:height)-1});
    % Plot occupancy map
    figure
    hold on
    hist3(dvtObjRun(:,3:4),'Edges',{(-width:gridSize:width)-1,(-height:gridSize:height)-1},'CDataMode','auto')
    scatter(objPosObjRun([3 5 7]), objPosObjRun([4 6 8]), 'y', 'filled')
    ylabel('Door Side')
    colorbar
    title('Object-centered')
    view(2)
    hold off

    % mapWorldCentered = occupancyMap(width,height);
    % mapObjCentered = occupancyMap(2*width,2*height);
    % 
    % [X,Y] = meshgrid(1:width, 1:height);
    % setOccupancy(mapWorldCentered, [X(:) Y(:)], hWorld(:))
    % setOccupancy(mapObjCentered, [X(:) Y(:)], hObj(:))
end

