%% Object locator
% Charles Xu @ UCSD, v1, 20221112
% Process scored csv file to get object location and events

%% Import scored csv file

[objFileName, objPathName] = uigetfile('*.csv', 'Choose the object marker file.');
[evtFileName, evtPathName] = uigetfile('*.csv', 'Choose the event marker file.');

objRaw = readtable(fullfile(objPathName,objFileName),'readvariablenames',false);
evtRaw = readtable(fullfile(evtPathName,evtFileName),'readvariablenames',false);

%% Save results - objVecMarkers

objVecMarkers.objPosition = objPosition;
objVecMarkers.

args = input('Save data? yes/no (y/n)','s');
if (args == "yes") | (args == 'y') %#ok<OR2>
    save(fullfile(dvtPathName,strcat(dvtFileName(1:end-4),'_RecStruct_ProcessedDVT_ObjVec')), 'indRecStruct');
end