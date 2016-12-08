fid = fopen('human_orthologs_top500scBATonly.csv');
top500_humanGenbank = {};
for i=1:2
    line = fgetl(fid);
end
j=1;

while line ~= -1 
    %line=strrep(line,'"','');
    words = strsplit(line,',');
    top500_humanGenbank(j,:)=[words(1),words(2)];
    words = {};
    line = {};
    line = fgetl(fid);
    j=j+1;
end
fclose(fid);

fid = fopen('human_orthologs_top500scBAT_2fold.csv');
top500_humanGenbank_2fold = {};
for i=1:2
    line = fgetl(fid);
end
j=1;

while line ~= -1 
    %line=strrep(line,'"','');
    words = strsplit(line,',');
    top500_humanGenbank_2fold(j,:)=[words(1),words(2)];
    words = {};
    line = {};
    line = fgetl(fid);
    j=j+1;
end
fclose(fid);
