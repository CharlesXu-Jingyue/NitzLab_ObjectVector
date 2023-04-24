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
close all

% Prompt user to select file. Will save back to the same folder
currectDir = pwd;
recDir = uigetdir;
cd(recDir)
[matFileName, matPathName] = uigetfile(fullfile(recDir, '*indRecStruct.mat'), 'Choose the mat file.');
load(fullfile(matPathName, matFileName))
figSaveDir = fullfile(string(recDir) + filesep + "PosTuning");
if ~exist(figSaveDir, 'dir')
   mkdir(figSaveDir)
end
cd(currectDir)

% Resolution
width = 640;
height = 480;

% Sampling rate
Fs = 60;

% Bin size in pixels
binSize = 20;

%% Analysis
iInner = indRecStruct.event.iInner;
tInner = indRecStruct.event.tInner;
frameRef = ["World", "Object"];
for i = 1:length(frameRef)
    % Read dvt and object position data
    if frameRef(i) == "World"
        dvt = indRecStruct.world.processedDVT(:,[1,2,9,10]);
        pos = dvt(:,3:4);
        pos(pos==1) = NaN;
        dvt(:,3:4) = pos;
        objPos = indRecStruct.world.objPosition;
    elseif frameRef(i) == "Object"
        dvt = indRecStruct.object.processedDVT(:,[1,2,9,10]);
        pos = dvt(:,3:4);
        pos(pos==1) = NaN;
        dvt(:,3:4) = pos;
        objPos = indRecStruct.object.objPosition;
    end
    
    % Read spiking data
    fieldNames = string(fieldnames(indRecStruct.spike));
    fieldCell = struct2cell(indRecStruct.spike);
    isSignalCell = strfind(fieldNames, "sig");
    isSignalCell(cellfun(@isempty, isSignalCell)) = {0};
    isSignal = logical([isSignalCell{:}]);
    signalNames = fieldNames(isSignal);
    variableCell = fieldCell(isSignal);
    
    numSpikes = cell(numel(variableCell),1);
    posTuning = cell(numel(variableCell),1);
    posInfo = zeros([numel(variableCell),2]);
    
    % Only include inner occupancy
    dvtInner = [];
    for j = 1:size(iInner)
        dvtInner = [dvtInner; dvt(iInner(j,2):iInner(j,3),:)]; 
    end
    
    if frameRef(i) == "World"
        occupancy = hist3(dvtInner(:,3:4), 'Edges', {0:binSize:width 0:binSize:height}, 'CDataMode', 'auto');
    elseif frameRef(i) == "Object"
        occupancy = hist3(dvtInner(:,3:4), 'Edges', {-width:binSize:width -height:binSize:height}, 'CDataMode', 'auto');
    end
    
    % tic
    for j = 1:size(variableCell, 1)    
        % Compute rate map for each cell
        jCell = variableCell{j,1};
        isInRange = (jCell >= tInner(:,2)') & (jCell <= tInner(:,3)');
        jCellInner = jCell(any(isInRange, 2));
        
        % Use sort
        sortCell = [(1:size(dvtInner,1))' dvtInner(:,2); 
            repelem(size(dvtInner,1)+numel(jCellInner)+1, numel(jCellInner))' jCellInner];
        sortCell = sortrows(sortCell, 2);
        indCellInner = (sortCell(:,1) - (1:size(sortCell, 1))') > 0;
        indCellInner = (nonzeros(indCellInner .* (1:size(sortCell, 1))') - (1:size(jCellInner, 1))');
        indCellAll = dvtInner(indCellInner,1);
        jCellInner = [indCellAll jCellInner]; %#ok<*AGROW>
        
        % Use min
%         for iSpk = 1:length(jCell)
%             hdDist = abs(jCell(iSpk,1)-dvtInner(:,2));
%             [~, nearestSpkInd] = min(hdDist);
%             SpkHd = dvtInner(nearestSpkInd,end);
%             jCell(iSpk,2) = SpkHd;
%         end
        
        % Spike per bin
        if frameRef(i) == "World"
            jOcc = hist3(dvtInner(indCellInner, 3:4), 'Edges', {0:binSize:width 0:binSize:height}, 'CDataMode', 'auto');
        elseif frameRef(i) == "Object"
            jOcc = hist3(dvtInner(indCellInner, 3:4), 'Edges', {-width:binSize:width -height:binSize:height}, 'CDataMode', 'auto');
        end
        numSpikes{j,1} = jOcc;
        
        % Positional tuning
        jPosTuning = jOcc ./ occupancy * Fs;
%         imagesc(flip(jPosTuning',1)); colorbar; grid on; xticks(0.5:1:size(occupancy,1)); yticks(0.5:1:size(occupancy,2));
        posTuning{j,1} = jPosTuning;

        % Calculate positional information per cell
        rm = jPosTuning;
        occ = occupancy;
        occThresh = 1;
        Ro = numel(jCell)/dvt(end,2); % Overall firing rate
        [infoPerSecond, infoPerSpike] = Doug_spatialInfo(rm,Ro,occ,occThresh);
        posInfo(j,:) = [infoPerSecond, infoPerSpike];
        
        % Plot spike counts, occupancy, and rate maps
        figure(1)
        clf
        set(gcf,'Position',[0   319   1682   294])
        subplot(1,3,1)
        imagesc(rot90(jOcc))
        colorbar
        title('Number of spikes')
        subplot(1,3,2)
        imagesc(rot90(occupancy))
        colorbar
        title('Occupancy')
        subplot(1,3,3)
        imagesc(rot90(jPosTuning))
        colorbar
        title('Tuning Curve (Hz)')

        saveas(gcf, figSaveDir + filesep + "Signal" + j + "_" + frameRef(i) + ".png")
        pause(1)
    end
    
    % Save occupancy matrix, save rate map, and save infoPerSecond and infoPerSpike for each cell to each row
    if frameRef(i) == "World"
        indRecStruct.world.posNumSpikes = numSpikes;
        indRecStruct.world.posOccupancy = occupancy;
        indRecStruct.world.posTuning = posTuning;
        indRecStruct.world.posInfo = posInfo;
    elseif frameRef(i) == "Object"
        indRecStruct.object.posNumSpikes = numSpikes;
        indRecStruct.object.posOccupancy = occupancy;
        indRecStruct.object.posTuning = posTuning;
        indRecStruct.object.posInfo = posInfo;
    end
    
    % toc
end

%% Save data
args = input('Save data? yes/no (y/n)','s');
if (args == "yes") | (args == 'y') %#ok<OR2>
    save(fullfile(matPathName, matFileName), 'indRecStruct');
end

