clear

recDir = uigetdir;
cd(recDir)
[matFileName, matPathName] = uigetfile(fullfile(recDir, '*indRecStruct*.mat'), 'Choose the mat file.');
load(fullfile(matPathName, matFileName))
fieldNames = string(fieldnames(indRecStruct.spike));
fieldCell = struct2cell(indRecStruct.spike);
isSignalCell = strfind(fieldNames, "sig");
isSignalCell(cellfun(@isempty, isSignalCell)) = {0};
isSignal = logical([isSignalCell{:}]);
signalNames = fieldNames(isSignal);
variableCell = fieldCell(isSignal);

save allTfiles.mat variableCell

fid = fopen('tfile_list.txt', 'w+');
for i = 1:length(signalNames)
    filename = signalNames(i,1);
    signal = variableCell{i};
    codethis = ['save ' sprintf(char(filename)) ' signal -ascii'];
    eval(codethis);
    fprintf(fid, '%s\n', char(filename));
end
fclose(fid);
