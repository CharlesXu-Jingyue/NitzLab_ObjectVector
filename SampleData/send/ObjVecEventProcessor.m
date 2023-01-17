%% ObjVecEventProcessor
% Charles Xu @ UCSD, v1, 20221114
% ObjVecEventProcessor reads the preprocessed indRecStruct data and event
% markers to process behavioral variables and task phases

%% Load processed indRecStruct file and event file
clear

[matFileName, matPathName] = uigetfile('*.mat', 'Choose the indRecStruct file.');
load(fullfile(matPathName,matFileName))
[evtFileName, evtPathName] = uigetfile('*.csv', 'Choose the event marker file.');
evtRaw = readtable(fullfile(evtPathName,evtFileName),'readvariablenames',false);

%% Process event markers
inner = table2array(evtRaw(evtRaw.(1)=="inner",5:6));
outer = zeros(length(inner)-1, 2);
outer(:,1) = inner(1:end-1, 2);
outer(:,2) = inner(2:end, 1);

indRecStruct.event.inner = inner;
indRecStruct.event.outer = outer;

%% Save results
args = input('Save data? yes/no (y/n)','s');
if (args == "yes") | (args == 'y') %#ok<OR2>
    save(fullfile(matPathName,strcat(matFileName(1:end-4),'_eventProcessed')), 'indRecStruct');
end