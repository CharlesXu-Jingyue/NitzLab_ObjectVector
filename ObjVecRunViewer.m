%% ObjVecRunViewer
% ObjVecRunViewer reads in the indRecStruct dataset and plots single runs
% 
% By Jingyue Xu, 202203317

clear
close all

% Prompt user to select file
recDir = uigetdir;
cd(recDir)
[matFileName, matPathName] = uigetfile(fullfile(recDir, '*indRecStruct.mat'), 'Choose the mat file.');
load(fullfile(matPathName, matFileName))

% Prompt user to specify run
%args = input("This recording has " + indRecStruct.event.runNumber + " runs, which run to plot?");
%
%iRun = args;
iRun = 15;
runIndex = indRecStruct.event.inner(iRun,2:3);

dvtWorld = indRecStruct.world.processedDVT(runIndex(1):runIndex(2),9:10);
dvtObject = indRecStruct.object.processedDVT(runIndex(1):runIndex(2),9:10);
objPosWorld = indRecStruct.world.objPosition(runIndex(1),3:end);
objPosObject = indRecStruct.object.objPosition(runIndex(1),3:end);

HDRWorld = indRecStruct.world.HDRadians(runIndex(1):runIndex(2),1);
HDRObject = indRecStruct.object.HDRadians(runIndex(1):runIndex(2),1);

distWorld = indRecStruct.world.objVec(runIndex(1):runIndex(2),4);
distObject = indRecStruct.object.objVec(runIndex(1):runIndex(2),4);

c = linspace(runIndex(1), runIndex(2), runIndex(2)-runIndex(1)+1);
runIndices = runIndex(1):1:runIndex(2);

% Plot raw trace
figure(1)
hold on
scatter(dvtWorld(:,1), dvtWorld(:,2), [], c)
scatter(objPosWorld(1,1:2:end), objPosWorld(1,2:2:end), 'r', 'filled')
title('Raw Trace in World Frame')
colorbar
hold off

figure(2)
hold on
scatter(dvtObject(:,1), dvtObject(:,2), [], c)
scatter(objPosObject(1,1:2:end), objPosObject(1,2:2:end), 'r', 'filled')
title('Raw Trace in Object Frame')
colorbar
hold off

% Plot HDR over trial
figure(3)
hold on
plot(runIndices, HDRWorld)
ylim([-pi, pi])
title('HDR over Run in World Frame')
hold off

figure(4)
hold on
plot(runIndices, HDRObject)
ylim([-pi, pi])
title('HDR over Run in Object Frame')
hold off

% Distance over time
figure(5)
hold on
plot(runIndices, distWorld)
title('Distance to Object')
hold off


