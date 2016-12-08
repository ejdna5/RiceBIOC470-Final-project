function [data] = readCSV_ChenData(filename)

fid = fopen(filename);
data = {};
line = fgetl(fid);
i=1;

while line ~= -1 
    words = strsplit(line,',');
    data(i,:)=[words(11),words(12),words(2),words(3),words(4),words(6)];
    words = {};
    line = {};
    line = fgetl(fid);
    i=i+1;
end
fclose(fid);
end

