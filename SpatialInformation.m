%% Selected cells
load('ratemap_struct_22.mat');

cell_nos = [54, 55, 59, 61, 62, 70, 88, 107, 115, 136];
for i=1:length(cell_nos)
    j=cell_nos(i);
    j
    
    rm = ratemap_struct(j).rm_ns; 
    Ro = nanmean(nanmean(rm));
    occ = ratemap_struct(j).occMap;
    occThresh = 1;
    [infoPerSecond, infoPerSpike] = Doug_spatialInfo(rm,Ro,occ,occThresh)
end
    

%% All of them
