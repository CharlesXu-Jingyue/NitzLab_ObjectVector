[csvFileName, csvPathName] = uigetfile('*.csv', 'Choose the marker file.');
rawMarker = readtable(fullfile(csvPathName,csvFileName),'readvariablenames',false);

a = readmatrix(fullfile(csvPathName,csvFileName),'FileType','text','OutputType','double');

rawMarker(:,1)
size(rawMarker,1)

for i = 1:size(rawMarker,1)
    if rawMarker(:,1) == "Lego"
        
    end
end

a = 'Lego';
class(a)

rem(10, 4)

rawMarker(rawMarker.(1)=="Lego",:)
