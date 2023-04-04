%% Selected Cells
load('ratemap_struct_22.mat');
cell_nos = [54, 55, 59, 61, 62, 70, 88, 107, 115, 136];

half_narrow_val = 1;
narrow_stdev_val = 0.5;
half_narrow = half_narrow_val;
narrow_stdev = narrow_stdev_val;
   
for i=1:length(cell_nos)
    j=cell_nos(i);
    
    ratemap = ratemap_struct(j).rm_ns;
    
    [x y]= meshgrid(-half_narrow:1:half_narrow);
    narrow_gaussian = double_gaussian(x,y,narrow_stdev);
    narrow_gaussian = narrow_gaussian./sum(sum(narrow_gaussian));
    ratemap_smoothed = nanconv(ratemap, narrow_gaussian);
    % ratemap_smoothed = conv2(squeeze(ratemap), narrow_gaussian);
    
    f = size(ratemap_smoothed);
    corr_cols = zeros(f(1)*f(2)*8,2);
    count = 0;
    for a = 2:f(1)-1
        for b = 2:f(2)-1
            corr_cols(count+1:count+8,1) = ratemap_smoothed(a,b);
            corr_cols(count+1,2) = ratemap_smoothed(a-1,b);
            corr_cols(count+2,2) = ratemap_smoothed(a+1,b);
            corr_cols(count+3,2) = ratemap_smoothed(a,b-1);
            corr_cols(count+4,2) = ratemap_smoothed(a,b+1);
            corr_cols(count+5,2) = ratemap_smoothed(a+1,b+1);
            corr_cols(count+6,2) = ratemap_smoothed(a-1,b-1);
            corr_cols(count+7,2) = ratemap_smoothed(a-1,b+1);
            corr_cols(count+8,2) = ratemap_smoothed(a+1,b-1);
            count=count+8;
        end
    end
    
    h=isnan(corr_cols(:,1));
    hh=isnan(corr_cols(:,2));
    h(:,2)=hh;
    hhh = zeros(length(h),1);
            
    for c=1:length(h)
        p = h(c,:);
        if sum(p)~=0
            hhh(c)=1;
        end
    end
    
    hhh = logical(hhh);
    corr_cols(hhh,:) = [];

    j
    z = corrcoef(corr_cols);
    coherence = z(1,2)
    figure(1)
    plot(corr_cols(:,1),corr_cols(:,2),'*')
    
    figure(2)
    y=ratemap_struct(j).rm_ns;
    imagesc(y)
    colorbar
    
    figure(3)
    imagesc(ratemap_smoothed)
    colorbar
    
    rm = ratemap; 
    Ro = nanmean(nanmean(rm));
    occ = ratemap_struct(j).occMap;
    occThresh = 1;
    [infoPerSecond, infoPerSpike] = Doug_spatialInfo(rm,Ro,occ,occThresh)

    pause
end

%% All Cells
clear
close all

load('ratemap_struct_22.mat');
coherence_vector = zeros(202, 1);
coherence_vector_ps = zeros(47, 1);
coherence_vector_ca = zeros(139, 1);
coherence_vector_su = zeros(16, 1);

spatial_info_vector = zeros(202, 1);
spatial_info_vector_ps = zeros(47, 1);
spatial_info_vector_ca = zeros(139, 1);
spatial_info_vector_su = zeros(16, 1);

ps_indx = 1;
ca_indx = 1;
su_indx = 1;

half_narrow_val = 1;
narrow_stdev_val = 0.5;
half_narrow = half_narrow_val;
narrow_stdev = narrow_stdev_val;

for i=1:202
    ratemap = ratemap_struct(i).rm_ns;
    cellRegion = ratemap_struct(i).cellRegion;
    
    [x y]= meshgrid(-half_narrow:1:half_narrow);
    narrow_gaussian = double_gaussian(x,y,narrow_stdev);
    narrow_gaussian = narrow_gaussian./sum(sum(narrow_gaussian));
    ratemap_smoothed = nanconv(ratemap, narrow_gaussian);
    % ratemap_smoothed = conv2(squeeze(ratemap), narrow_gaussian);
    
    f = size(ratemap_smoothed);
    corr_cols = zeros(f(1)*f(2)*8,2);
    count = 0;
    for a = 2:f(1)-1
        for b = 2:f(2)-1
            corr_cols(count+1:count+8,1) = ratemap_smoothed(a,b);
            corr_cols(count+1,2) = ratemap_smoothed(a-1,b);
            corr_cols(count+2,2) = ratemap_smoothed(a+1,b);
            corr_cols(count+3,2) = ratemap_smoothed(a,b-1);
            corr_cols(count+4,2) = ratemap_smoothed(a,b+1);
            corr_cols(count+5,2) = ratemap_smoothed(a+1,b+1);
            corr_cols(count+6,2) = ratemap_smoothed(a-1,b-1);
            corr_cols(count+7,2) = ratemap_smoothed(a-1,b+1);
            corr_cols(count+8,2) = ratemap_smoothed(a+1,b-1);
            count=count+8;
   
        end
    end
    
    h=isnan(corr_cols(:,1));
    hh=isnan(corr_cols(:,2));
    h(:,2)=hh;
    hhh = zeros(length(h),1);
            
    for c=1:length(h)
        p = h(c,:);
        if sum(p)~=0
            hhh(c)=1;
        end
    end
    
    hhh = logical(hhh);
    corr_cols(hhh,:) = [];

    %i
    z = corrcoef(corr_cols);
    coherence = z(1,2);
    
    coherence_vector(i) = coherence; 
    if cellRegion == 'ps'
        coherence_vector_ps(ps_indx) = coherence;
    end
    
    if cellRegion == 'ca'
        coherence_vector_ca(ca_indx) = coherence;
    end
    
    if cellRegion == 'su'
        coherence_vector_su(su_indx) = coherence;
    end
    
    
%     figure(1)
%     plot(corr_cols(:,1),corr_cols(:,2),'*')
%     
%     figure(2)
%     y=ratemap_struct(i).rm_ns;
%     imagesc(y)
%     colorbar
%     
%     figure(3)
%     imagesc(ratemap_smoothed)
%     colorbar
    
    rm = ratemap;
    Ro = nanmean(nanmean(rm));
    occ = ratemap_struct(i).occMap;
    occThresh = 1;
    [infoPerSecond, infoPerSpike] = Doug_spatialInfo(rm,Ro,occ,occThresh);
    
    spatial_info_vector(i) = infoPerSpike;
    if cellRegion == 'ps'
        spatial_info_vector_ps(ps_indx) = infoPerSpike;
        ps_indx = ps_indx + 1;
    end
    
    if cellRegion == 'ca'
        spatial_info_vector_ca(ca_indx) = infoPerSpike;
        ca_indx = ca_indx + 1;
    end
    
    if cellRegion == 'su'
        spatial_info_vector_su(su_indx) = infoPerSpike;
        su_indx = su_indx + 1;
    end
    
    
end

% Coherence vs. Spatial Information
figure()
plot(coherence_vector, spatial_info_vector, '*')
xlabel('Coherence')
ylabel('InfoPerSpike')
axis([0.1 0.8 0 6])
title('All regions')

figure()
% hist(spatial_info_vector_ca, 30)
plot(coherence_vector_ca, spatial_info_vector_ca, '*')
xlabel('Coherence')
ylabel('InfoPerSpike')
axis([0.1 0.8 0 6])
title('CA')

figure()
% hist(spatial_info_vector_su, 30)
plot(coherence_vector_su, spatial_info_vector_su, '*')
xlabel('Coherence')
ylabel('InfoPerSpike')
axis([0.1 0.8 0 6])
title('Subiculum')

figure()
plot(coherence_vector_ps, spatial_info_vector_ps, '*')
% hist(spatial_info_vector_ps, 30)
xlabel('Coherence')
ylabel('InfoPerSpike')
axis([0.1 0.8 0 6])
title('Parietal')
