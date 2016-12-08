fid = fopen('Shinoda.txt');
humanScBATvsWAT = {};
for i=1:8
    line = fgetl(fid);
end
j=1;

while line ~= -1 
    line=strrep(line,'"','');
    words = strsplit(line,' ');
    humanScBATvsWAT(j,:)=[words(1),words(2),words(end-2),words(end-1),words(end)];
    words = {};
    line = {};
    line = fgetl(fid);
    j=j+1;
end
fclose(fid);

