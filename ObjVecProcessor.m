%% ObjVecProcessor
% ObjVecProcessor takes in a dvt file, the event csv files, and a spike
% data mat file and writes a .mat file with a name starting with the name
% of the dvt file.
% 
% The .mat file has interpolated light data up to acceptable gap sizes as
% well as an added "light" that is the average of the first two lights.
%
% This also includes velocity, acceleration, and head direction.
%
% Non-built-in functions called:
%      inpaint_nans
%
% By Jingyue Xu, 202203317
% Adapted from TTTrackingPreprocessor.m by Jake Olson, October 2014

clear

% Prompt user to select file. Will save back to the same folder.
recDir = uigetdir;
cd(recDir)
[dvtFileName, dvtPathName] = uigetfile(fullfile(recDir, '*.dvt'), 'Choose the dvt file.');
[objFileName, objPathName] = uigetfile(fullfile(recDir, '*lego.csv'), 'Choose the object marker file.');
[inrFileName, inrPathName] = uigetfile(fullfile(recDir, '*inner.csv'), 'Choose the inner marker file.');
[spkFileName, spkPathName] = uigetfile(fullfile(recDir, '*.mat'), 'Choose the spike data file.');
indRecStruct.dvtFileName = dvtFileName;
indRecStruct.dvtPathName = dvtPathName;

% Raw DVT file should be of the format of each row being a sample, with
% multiple columns. The columns are:
% Line Time Light1X Light1Y Light2X Light2Y ...
dvtRaw = load(fullfile(dvtPathName, dvtFileName));
workingDVT = dvtRaw;

% Read object location markers from CSV file
objRaw = readtable(fullfile(objPathName, objFileName));

% Read inner run event markers from CSV file
inrRaw = readtable(fullfile(inrPathName, inrFileName));

% Read spike time data MAT file
spkRaw = load(fullfile(spkPathName, spkFileName));

% So we can use the tracked pixel values as indices later on in analyses,
% we add one. To encode the location of the light in the DVT files, Plexon
% uses values from a range of 1-1024? w/ 0 values indicating lost tracking.
% This converts from the encoding form 1024x1024 to native pixels of camera
% 640x480, then adds 1 to the value so values are 1-640 instead of 0-639.
workingDVT(:, 3:end) = workingDVT(:, 3:end)*639/1023+1;
nRealLights = (size(dvtRaw, 2)/2)-1;

clear objFileName objPathName inrFileName inrPathName spkFileName spkPathName

%% Filling in missing position points.
% For this work, Doug is willing to fill in gaps of up to 1/2 second.
% I use points 99 and 100 just in case something strange happens in the
% startup. Sufficiently large to avoid that IMO.
sampleRate = round(1/(workingDVT(100, 2)-workingDVT(99, 2))); 
maxGap = sampleRate/2; % Half second worth of samples.

indRecStruct.trackingSampleRate = sampleRate;
indRecStruct.maxGapFilled = maxGap;

% Var Init
samplesLost = false(length(dvtRaw), nRealLights);
samplesFilled = false(length(dvtRaw), nRealLights);
samplesUnfilled = false(length(dvtRaw), nRealLights);

% Run for each light.
for iLight = 1:nRealLights
    
    %Indices of the light columns in the dvt matrix.
    lightColX = 1+iLight*2;
    lightColY = 2+iLight*2;
    
    % Find the edges of the gaps.
    lostTrackingEdges = diff([0;workingDVT(:, lightColX) == 1] &...
        [0;workingDVT(:, lightColY) == 1]);
    
    % Mark the bad points w/ NaNs.
    workingDVT(...
        (workingDVT(:,lightColX) == 1 & workingDVT(:,lightColY) == 1),...
        lightColX:lightColY) = NaN;
    samplesLost(:,iLight) = isnan(workingDVT(:,lightColX));
    
    % Interpolate values - using the inpaint_nans fn. Rounded because the
    % values are indices of pixels. Our max precision is individual pixels.
    workingDVT(:,lightColX:lightColY) = ...
        round(inpaint_nans(workingDVT(:,lightColX:lightColY)));
    
    % Find the gaps too large to fix.
    gapStarts = find(lostTrackingEdges == 1);
    gapEnds = find(lostTrackingEdges == -1);
    if ~isempty(gapStarts)
%         if gapEnds(1) < gapStarts(1) % Lost tracking at beginning. Common.
%             gapStarts = [1;gapStarts];
%         end
        if gapStarts(end) > gapEnds(end) % Lost tracking at end. Common.
            gapEnds = [gapEnds;length(workingDVT)+1]; %#ok<AGROW>
        end
        gapLengths = gapEnds - gapStarts;
        
        % "Unfix" the gaps longer than the length we are willing to
        % interpolate over - Done this way so that when the interpolation
        % algorithm runs, there are no bad data (minus spurious
        % reflections, etc.) that could bias the algorithm in a certain
        % direction.
        unfixableGaps = gapLengths > maxGap;
        for iGap = 1:length(gapStarts)
            if unfixableGaps(iGap)
                workingDVT(gapStarts(iGap):gapEnds(iGap)-1,lightColX:lightColY) = 1;
                % [1,1] light coordinate now represents lost tracking.
            end
        end
    end
    
    % So we know where we fixed.
    samplesFilled(:,iLight) = samplesLost(:,iLight) & ...
        (workingDVT(:,lightColX) ~= 1 | workingDVT(:,lightColY) ~= 1);
    
    % The gaps we couldn't fix.
    samplesUnfilled(:,iLight) = samplesLost(:,iLight) & ...
        (workingDVT(:,lightColX) == 1 & workingDVT(:,lightColY) == 1);
end

indRecStruct.samplesLost = samplesLost;
indRecStruct.samplesFilled = samplesFilled;
indRecStruct.samplesUnfilled = samplesUnfilled;

clear gap* iGap iLight lightCol* lostTrackingEdges unfixableGaps maxGap

%% Create a "light" and add to DVT matrix that is the average of the first two lights.
% Var Init
avgLight = nRealLights+1; % Adding the average of the lights as a 'light'.

% Average the XY coords of the lights for the entire recording - save as 2
% new columns - a new *light* for the DVT matrix.
workingDVT(:,(1+avgLight*2)) = round(sum(workingDVT(:,3:2:1+nRealLights*2),2)/nRealLights);
workingDVT(:,(2+avgLight*2)) = round(sum(workingDVT(:,4:2:2+nRealLights*2),2)/nRealLights);

% Put [1,1] (missing light code) into the averaged columns for
% samples where one or more of the lights is lost.
notPerfectTracking = any(samplesUnfilled,2);
workingDVT(notPerfectTracking,(1+avgLight*2):(2+avgLight*2)) = 1;

clear iLight lightCol* avgLight

%% Create a mashup "light" 
% Add to DVT matrix a "light" that is the average of the
% first two lights or just the value of each individual light if the other 
% is lost.

% Var Init
mashupLight = nRealLights+2; % Adding the average of the lights as a 'light'.

% Average the XY coords of the lights for the entire recording - save as 2
% new columns - a new *light* for the DVT matrix.
workingDVT(:,(1+mashupLight*2)) = round(sum(workingDVT(:,3:2:1+nRealLights*2),2)/nRealLights);
workingDVT(:,(2+mashupLight*2)) = round(sum(workingDVT(:,4:2:2+nRealLights*2),2)/nRealLights);

% Put [1,1] (missing light code) into the averaged columns for
% samples where one or more of the lights is lost.
workingDVT(any(samplesUnfilled,2),(1+mashupLight*2):(2+mashupLight*2)) = 1;

% Fill spots where we can't average but do have 1 light (lost 1 light).
for iLight = 1:nRealLights
    % Indices of the light columns in the dvt matrix.
    lightColX = 1+iLight*2;
    lightColY = 2+iLight*2;
    thisLightIsGood = ~samplesUnfilled(:,iLight);
    workingDVT(thisLightIsGood & notPerfectTracking,(1+mashupLight*2):(2+mashupLight*2)) =...
        workingDVT(thisLightIsGood & notPerfectTracking,lightColX:lightColY);
end

% DVT processing - (interpolating and adding an average light) is finished.
processedDVT = workingDVT;
indRecStruct.world.processedDVT = processedDVT;

clear iLight lightCol* workingDVT mashupLight thisLightIsGood notPerfectTracking

%% Process object markers into object location data
objPosition = zeros(size(processedDVT,1),6); % Getting object positions throughout recording

rawObjMarker = objRaw(objRaw.(1) == "lego" | objRaw.(1) == "Lego" ,:);
for i = (1:size(rawObjMarker,1)/3)-1
    tMarker = rawObjMarker{3*i+1,5};
    objPosition1 = rawObjMarker{3*i+1,7:8};
    objPosition2 = rawObjMarker{3*i+2,7:8};
    objPosition3 = rawObjMarker{3*i+3,7:8};
    iMarker = round(tMarker*60);
    objPosition(iMarker:end,1) = objPosition1(1)*640/1024; % columns 1, 2 are x, y positions for the light at arm A (right side of the angle) of the object
    objPosition(iMarker:end,2) = objPosition1(2)*480/768;
    objPosition(iMarker:end,3) = objPosition2(1)*640/1024; % columns 3, 4 are x, y positions for the light at the vertex of the object
    objPosition(iMarker:end,4) = objPosition2(2)*480/768;
    objPosition(iMarker:end,5) = objPosition3(1)*640/1024; % columns 5, 6 are x, y positions for the light at arm B (left side of the angle) of the object
    objPosition(iMarker:end,6) = objPosition3(2)*480/768;
end

objPosition = [processedDVT(:,1:2) objPosition];
indRecStruct.world.objPosition = objPosition;

clear tMarker objPosition1 objPosition2 objPosition3 iMarker

%% Create relative DVT for object-centered position
% Apply a transformation matrix to rotate DVT to object-relative position
% This code block is written assuming there are 3 lights (A, B, C) on the
% object, where B is at the corner of the object, A is on the right to B

workingDVTRel = processedDVT;
objPositionRel = objPosition(:,3:end);

for i = 1:size(workingDVTRel,2)/2-1
    workingDVTRel(:,i*2+1:i*2+2) = workingDVTRel(:,i*2+1:i*2+2) - objPosition(:,5:6);
end

for i = 1:size(objPositionRel,2)/2
    objPositionRel(:,i*2-1:i*2) = objPositionRel(:,i*2-1:i*2) - objPosition(:,5:6);
end

for i = 1:size(workingDVTRel,1)
    vecAxy = objPositionRel(i,1:2) - objPositionRel(i,3:4); % Vector A corresponding to object x-axis
    theta = -atan2(vecAxy(2),vecAxy(1)); % Object angle relative to room coordinates
    A = [cos(theta), -sin(theta); sin(theta), cos(theta)]; % Create a rotation matrix
    
    for j = 1:size(workingDVTRel,2)/2-1
        if workingDVTRel(i,j*2+1) ~= 1 && workingDVTRel(i,j*2+1) ~= 1
            workingDVTRel(i,j*2+1:j*2+2) = (A*workingDVTRel(i,j*2+1:j*2+2)')'; % Convert into object-relative coordinates 
        end
    end
    
    for j = 1:size(objPositionRel,2)/2
        objPositionRel(i,j*2-1:j*2) = (A*objPositionRel(i,j*2-1:j*2)')';
    end
end

% for i = 1:size(workingDVTRel,2)/2-1
%     workingDVTRel(:,i*2+1:i*2+2) = workingDVTRel(:,i*2+1:i*2+2) + objPosition(:,3:4);
% end
% 
% for i = 1:size(objPositionRel,2)/2
%     objPositionRel(:,i*2-1:i*2) = objPositionRel(:,i*2-1:i*2) + objPosition(:,3:4);
% end

objPositionRel = [processedDVT(:,1:2) objPositionRel];
indRecStruct.object.processedDVT = workingDVTRel;
indRecStruct.object.objPosition = objPositionRel;

clear vecAxy theta

%% Calculate absolute distance to object (vertex)
distance = zeros(size(processedDVT,1),4);
distance(:,1:2) = processedDVT(:,1:2);
distance(:,3) = sqrt((processedDVT(:,9)-objPosition(:,3)).^2 + (processedDVT(:,10)-objPosition(:,4)).^2);
posVec = processedDVT(:,9:10)-objPosition(:,3:4);
distance(:,4) = atan2(posVec(:,2), posVec(:,1));

indRecStruct.world.distanceToObj = distance;

clear posVec

%% Velocity & Acceleration - Averaged over adaptable window (updated with object-centered direction)
% Initialize window size to use - can change here if desired.
velSmoothWinSecs = 1/10; % Uses position change over X sec to calc vel.
velSmoothWinSamples = round(sampleRate*velSmoothWinSecs);

% Output Variable Init
vel = nan(length(processedDVT)-velSmoothWinSamples,2,nRealLights);
acc = nan(length(vel)-1,2,nRealLights); % Compare point by point vel since they are already smoothed.
instVel = nan(length(processedDVT)-1,2,nRealLights); % Instantaneous (sample rate) velocity
instAcc = nan(length(instVel)-1,2,nRealLights); % Compare point by point vel since they are already smoothed.

velRel = vel;
accRel = acc;
instVelRel = instVel;
instAccRel = instAcc;

vecAxy = objPosition(:,1:2) - objPosition(:,3:4); % Get array of object directions
theta = -atan2(vecAxy(:,2),vecAxy(:,1));

for iLight = 1:nRealLights
    lightColX = 1+iLight*2;
    lightColY = 2+iLight*2;
    
    lightLow = processedDVT(1:end-velSmoothWinSamples,lightColX:lightColY);
    lightLow(repmat(samplesUnfilled(1:end-velSmoothWinSamples,iLight),1,2)) = NaN;
    lightHigh = processedDVT(1+velSmoothWinSamples:end,lightColX:lightColY);
    lightHigh(repmat(samplesUnfilled(1+velSmoothWinSamples:end,iLight),1,2)) = NaN;
    xyVel = lightHigh - lightLow;
    
    velMag = sqrt(sum((xyVel.^2),2))/velSmoothWinSecs; % In pixels/second.
    velDirection = atan2(xyVel(:,2),xyVel(:,1)); % -pi : pi
    vel(:,:,iLight) = [velMag,velDirection];
    xyAcc = diff(xyVel);
    
    accMag = sqrt(sum((xyAcc.^2),2))*sampleRate; % In pixels/second^2.
    accDirection = atan2(xyAcc(:,2),xyAcc(:,1)); % -pi : pi
    acc(:,:,iLight) = [accMag,accDirection];
    
    instXYVel = diff(processedDVT(:,lightColX:lightColY));
    instVelMag = sqrt(sum((instXYVel.^2),2))*sampleRate; % In pixels/second.
    instVelDirection = atan2(instXYVel(:,2),instXYVel(:,1)); % -pi : pi
    instVel(:,:,iLight) = [instVelMag,instVelDirection];
    instXYAcc = diff(instXYVel);
    instAccMag = sqrt(sum((instXYAcc.^2),2))*sampleRate; % In pixels/second.
    instAccDirection =  atan2(instXYAcc(:,2),instXYAcc(:,1)); % -pi : pi
    instAcc(:,:,iLight) = [instAccMag,instAccDirection];
    
    % Store relative orientations
    velDirectionRel = velDirection + theta(1:length(velDirection));
    velRel(:,:,iLight) = [velMag,velDirectionRel];
    accDirectionRel = accDirection + theta(1:length(accDirection));
    accRel(:,:,iLight) = [accMag,accDirectionRel];
    instVelDirectionRel = instVelDirection + theta(1:length(instVelDirection));
    instVelRel(:,:,iLight) = [instVelMag,instVelDirectionRel];
    instAccDirectionRel = instAccDirection + theta(1:length(instAccDirection));
    instAccRel(:,:,iLight) = [instAccMag,instAccDirectionRel];
end

indRecStruct.world.velInst = instVel;
indRecStruct.world.accInst = instAcc;
indRecStruct.world.velSmoothed = vel;
indRecStruct.world.accSmoothed = acc;

indRecStruct.object.velInst = instVelRel;
indRecStruct.object.accInst = instAccRel;
indRecStruct.object.velSmoothed = velRel;
indRecStruct.object.accSmoothed = accRel;

clear bufferDistance iPos iLight light* xyDiff speed*  velMag ...
    velDirection accMag accDirection

%% Head Direction
light1 = processedDVT(:,3:4);
light2 = processedDVT(:,5:6);

posDiff = light1-light2;
HDRadians = atan2(posDiff(:,2),posDiff(:,1));
HDRadians((samplesLost(:,1) & ~samplesFilled(:,1)) | ...
    (samplesLost(:,2) & ~samplesFilled(:,2))) = NaN;

indRecStruct.world.HDRadians = HDRadians;

% Relative Head Direction
HDRadiansRel = HDRadians + theta(1:length(HDRadians));
indRecStruct.object.HDRadians = HDRadiansRel;

% Output Check Code
% [count,center] = hist(HDRadians,36);
% sortedRows = sortrows([count;center]',1);

clear posDiff light* vecAxy theta

%% Process inner run event markers
inner = table2array(inrRaw(inrRaw.(1) == "inner" | inrRaw.(1) == "Inner", 5:6));
outer = zeros(length(inner)-1, 2);
outer(:,1) = inner(1:end-1, 2);
outer(:,2) = inner(2:end, 1);

indRecStruct.event.inner = inner;
indRecStruct.event.outer = outer;

clear inner outer

%% Process spike time data
indRecStruct.spike = spkRaw;

%% Save results - indRecStruct
args = input('Save data? yes/no (y/n)','s');
if (args == "yes") | (args == 'y') %#ok<OR2>
    save(fullfile(dvtPathName,strcat(dvtFileName(1:end-4),'_indRecStruct_objVec')), 'indRecStruct');
end

%% Visualize tracking and object position data
args = input('Plot data? yes/no (y/n)','s');
if (args == "yes") | (args == 'y') %#ok<OR2>
    figure
    hold on
    scatter(indRecStruct.world.processedDVT(:,9), indRecStruct.world.processedDVT(:,10), '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    scatter(indRecStruct.world.objPosition(1,[1,3,5]), indRecStruct.world.objPosition(1,[2,4,6]), 200, '.', 'MarkerEdgeColor', '#D95319')
    plot(indRecStruct.world.objPosition(1,[1,3,5]), indRecStruct.world.objPosition(1,[2,4,6]), 'LineWidth', 1, 'Color', '#D95319')
    hold off
    
    args = input('Plot object relative data? yes/no (y/n)','s');
    % scatter run and object data
    if (args == "yes") | (args == 'y') %#ok<OR2>
        figure
        hold on
        scatter(indRecStruct.object.processedDVT(:,9), indRecStruct.object.processedDVT(:,10), '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
        scatter(indRecStruct.object.objPosition(1,[1,3,5]), indRecStruct.object.objPosition(1,[2,4,6]), 200, '.', 'MarkerEdgeColor', '#D95319')
        plot(indRecStruct.object.objPosition(1,[1,3,5]), indRecStruct.object.objPosition(1,[2,4,6]), 'LineWidth', 1, 'Color', '#D95319')
        hold off
    end
end

%% Notes

% -[x] Make three lights on the object
% -[x] Shift numbers to positive
% -[x] Calculate distance to object
% -[x] Draw positional vectors with angle at some interval
% -[x] let the data tell you

% # Events
% ## Internal run
% Hard code gate locations
% Find stop events which correspond to getting the reward
% Find the previous crossing in and the next crossing out to define internal run
% 
% ## External run
% External run in between
% Clean/dirty run
% 
% # Maps
% ## Linearized rate map
% Need neuron data (20221114)
% 
% ## Rate map mapped against distance (need neuron data)
% At particular distance and orientation, how many spikes occurred
% Need neuron datva (20221114)
% 
% ### Occupancy map
% Matrix of zeros of the same size as the tracking, always centering at the object
% Update to get map for each inner run (20221114)
% 
% ### Spike map
% Based on occupancy map, how many spike occurred in a 1/60 s time window
% Need neuron data (20221114)


% Deal with cell phone
% Count unfixable gaps