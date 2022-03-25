function TTTrackingPreprocessor
%TTTrackingProcessor  2D Tracking PostProcessing
%   TTTRACKINGPROCESSOR takes in a dvt file and writes a .mat file with a
%   name starting with the name of the dvt file. The .mat file has
%   interpolated light data up to acceptable gap sizes as well as an
%   added "light" that is the average of the first two lights.
%
%   This also includes velocity, acceleration, and head direction.
%
%
%   Non-built-in functions called:
%        inpaint_nans
%
%   Written by Jake Olson, October 2014
%   Updated by Jake Olson, January 2015 - Vel, Acc, HD finished.
%   Last updated by Jake Olson, June 2016 - Changed output to struct.

% Prompt user to select file. Will save back to the same folder.
[dvtFileName, dvtPathName] = uigetfile('*.dvt', 'Choose the dvt file.');

indRecStruct.dvtFileName = dvtFileName;
indRecStruct.dvtPathName = dvtPathName;

% Raw DVT file should be of the format of each row being a sample, with
% multiple columns. The columns are:
% Line Time Light1X Light1Y Light2X Light2Y ...
rawDVT = load(fullfile(dvtPathName,dvtFileName));
workingDVT = rawDVT;

% So we can use the tracked pixel values as indices later on in analyses,
% we add one. To encode the location of the light in the DVT files, Plexon
% uses values from a range of 1-1024? w/ 0 values indicating lost tracking.
% This converts from the encoding form 1024x1024 to native pixels of camera
% 640x480, then adds 1 to the value so values are 1-640 instead of 0-639.
workingDVT(:,3:end) = workingDVT(:,3:end)*639/1023 +1;
nRealLights = (size(rawDVT,2)/2)-1;


%% Filling in missing position points.
% For this work, Doug is willing to fill in gaps of up to 1/2 second.
% I use points 99 and 100 just in case something strange happens in the
% startup. Sufficiently large to avoid that IMO.
sampleRate = round(1/(workingDVT(100,2)-workingDVT(99,2))); 
maxGap = sampleRate/2; % Half second worth of samples.

indRecStruct.trackingSampleRate = sampleRate;
indRecStruct.maxGapFilled = maxGap;

% Var Init
samplesLost = false(length(rawDVT),nRealLights);
samplesFilled = false(length(rawDVT),nRealLights);
samplesUnfilled = false(length(rawDVT),nRealLights);

% Run for each light.
for iLight = 1:nRealLights
    %Indices of the light columns in the dvt matrix.
    lightColX = 1+iLight*2;
    lightColY = 2+iLight*2;
    
    % Find the edges of the gaps.
    lostTrackingEdges = diff([0;workingDVT(:,lightColX) == 1] &...
        [0;workingDVT(:,lightColY) == 1]);
    
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
            gapEnds = [gapEnds;length(workingDVT)+1];
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

clear gap* iGap iLight lightCol* lostTrackingEdges unfixableGaps

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
    %Indices of the light columns in the dvt matrix.
    lightColX = 1+iLight*2;
    lightColY = 2+iLight*2;
    thisLightIsGood = ~samplesUnfilled(:,iLight);
    workingDVT(thisLightIsGood & notPerfectTracking,(1+mashupLight*2):(2+mashupLight*2)) =...
        workingDVT(thisLightIsGood & notPerfectTracking,lightColX:lightColY);
end

% DVT processing - (interpolating and adding an average light) is finished.
processedDVT = workingDVT;
indRecStruct.processedDVT = processedDVT;

clear iLight lightCol* workingDVT mashupLight thisLightIsGood notPerfectTracking

%% Velocity & Acceleration - Averaged over adaptable window.
% Initialize window size to use - can change here if desired.
velSmoothWinSecs = 1/10; %Uses position change over X sec to calc vel.
velSmoothWinSamples = round(sampleRate*velSmoothWinSecs);

indRecStruct.velSmoothWinSecs = velSmoothWinSecs;

% Output Variable Init
vel = nan(length(processedDVT)-velSmoothWinSamples,2,nRealLights);
acc = nan(length(vel)-1,2,nRealLights); %Compare point by point vel since they are already smoothed.

velInst = nan(length(processedDVT)-1,2,nRealLights); % instantaneous (sample rate) velocity
accInst = nan(length(velInst)-1,2,nRealLights); %Compare point by point vel since they are already smoothed.


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
end

indRecStruct.velInst = instVel;
indRecStruct.accInst = instAcc;
indRecStruct.velSmoothed = vel;
indRecStruct.accSmoothed = acc;

clear bufferDistance iPos iLight light* xyDiff speed* 
clear velDirection velMag accDirection accMag

%% Head Direction
light1 = processedDVT(:,3:4);
light2 = processedDVT(:,5:6);

posDiff = light1-light2;
HDRadians = atan2(posDiff(:,2),posDiff(:,1));
HDRadians((samplesLost(:,1) & ~samplesFilled(:,1)) | ...
    (samplesLost(:,2) & ~samplesFilled(:,2))) = NaN;

indRecStruct.HDRadians = HDRadians;

% Output Check Code
% [count,center] = hist(HDRadians,36);
% sortedRows = sortrows([count;center]',1);

clear posDiff light*

%% Save results - indRecStruct
save(fullfile(dvtPathName,strcat(dvtFileName(1:end-4),'_RecStruct_ProcessedDVT')), 'indRecStruct');

end

