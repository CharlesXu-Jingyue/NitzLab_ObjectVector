%% ObjVecPosTuning
% ObjVecPosTuning performs positional tuning analysis, including the
% occupancy map, rate map, and positional information.
%
% Input: indRecStruct (post-processing)
% Output: occupancy, rate map, and positional information for each cell,
% saved in the indRecStruct
%
% By Charles Xu @ UCSD, 20230409
%
%% Initialize
clear

% Prompt user to select file. Will save back to the same folder
recDir = uigetdir;
cd(recDir)
[matFileName, matPathName] = uigetfile(fullfile(recDir, '*indRecStruct.mat'), 'Choose the mat file.');
load(fullfile(matPathName, matFileName))
figSaveDir = fullfile(string(recDir) + filesep + "HDTuning");
if ~exist(figSaveDir, 'dir')
   mkdir(figSaveDir)
end