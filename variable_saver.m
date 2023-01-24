save sig;

variableCell = load('sig.mat');
signalNames = fieldnames(variableCell);
fid = fopen('tfile_list.txt', 'w+');
for i = 1:length(signalNames)
    filename = signalNames(i,1);
    codethis = ['save ' sprintf(char(filename)) ' ' sprintf(char(filename)) ' -ascii'];
    eval(codethis);
    fprintf(fid, '%s\n',char(filename));
end
fclose(fid);
clear variableCell signalNames filename codethis;