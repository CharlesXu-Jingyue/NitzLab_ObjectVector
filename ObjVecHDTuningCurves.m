%% ObjVecHDTuningCurves
% Reads in the indRecStruct produced by ObjVecProcessor to generate head
% direction tuning curves for each spiking neuron.
% 
% By Charles Xu, 20230221

%% Initialize
clear

% Prompt user to select file. Will save back to the same folder.
recDir = uigetdir;
cd(recDir)
[matFileName, matPathName] = uigetfile(fullfile(recDir, '*indRecStruct*.mat'), 'Choose the mat file.');
load(fullfile(matPathName, matFileName))

% Prompt user to specify frame of reference
args = input('Which frame of reference? (world/object)', 's');
if args == "world"
    processedDVT = indRecStruct.world.processedDVT;
    HDRadians = indRecStruct.world.HDRadians;
elseif args == "object"
    processedDVT = indRecStruct.object.processedDVT;
    HDRadians = indRecStruct.object.HDRadians;
else
    print('Please specify world or object')
end

% Read spiking data
fieldNames = string(fieldnames(indRecStruct.spike));
fieldCell = struct2cell(indRecStruct.spike);
isSignalCell = strfind(fieldNames, "sig");
isSignalCell(cellfun(@isempty, isSignalCell)) = {0};
isSignal = logical([isSignalCell{:}]);
signalNames = fieldNames(isSignal);
variableCell = fieldCell(isSignal);

% video-tracking sampling frequency
Fs = 60;

%angular bins
da = pi/30; %6 degrees
angBins = [-pi+da/2:da:pi-da/2];
ang = [processedDVT(:,2),HDRadians(:,1)];

%Occupancy
histAng = hist(ang(:,2),angBins);

% Which cell to plot
for cellIx = 1:length(variableCell)
    spk = variableCell{cellIx};
    
    for iSpk = 1:length(spk)
        
        hdDist = abs(spk(iSpk,1)-ang(:,1));
        [~, nearestSpkInd] = min(hdDist);
        SpkHd = ang(nearestSpkInd,end);
        spk(iSpk,2) = SpkHd;
        
    end
    
    
    % histogram of the number of times the cell fired in each bin of
    % head-direction
    spkPerAng = hist(spk(:,2),angBins);
    
    % now compute the tuning curve:
    hdTuning = spkPerAng./histAng * Fs;
    
    figure(1),clf
    set(gcf,'Position',[62   319   783   281])
    subplot(1,3,1)
    polarplot(angBins,spkPerAng)
    title('Number of spikes')
    subplot(1,3,2)
    polarplot(angBins,histAng)
    title('Occupancy')
    subplot(1,3,3)
    polarplot(angBins,hdTuning)
    title('Tuning Curve (Hz)')
    
    pause;
end
