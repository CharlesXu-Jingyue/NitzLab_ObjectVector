[csvFileName, csvPathName] = uigetfile('*.csv', 'Choose the marker file.');
rawMarker = readtable(fullfile(csvPathName,csvFileName),'readvariablenames',false);

