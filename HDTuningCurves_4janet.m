% Example script showing how to compute the head-direction tuning curve
% for single dataset


% File inputs:
% allTfiles: a cell array that contains 19 elements (19 HD cells
% simultaneously recorded) 
% processedDVT: sampleCount,timeStamp,LED-matrix
% HDRadians: 1 column matrix (ang -pi:pi radians)

%% Input setup

Fs = 60; % video-tracking sampling frequency

%angular bins
da = pi/30; %6 degrees
angBins = [-pi+da/2:da:pi-da/2];
ang = [processedDVT(:,2),HDRadians(:,1)];

%Occupancy
histAng = hist(ang(:,2),angBins);

% Which cell to plot
for cellIx = 1:length(allTfiles)
    spk = allTfiles{cellIx};
    
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
    polar(angBins,spkPerAng)
    title('Number of spikes')
    subplot(1,3,2)
    polar(angBins,histAng)
    title('Occupancy')
    subplot(1,3,3)
    polar(angBins,hdTuning)
    title('Tuning Curve (Hz)')
    
    pause;
    
end
