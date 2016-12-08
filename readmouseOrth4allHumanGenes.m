fid = fopen('mouse_orthologs_of_human_genes.csv');
humanGeneIDconversion = {};
for i=1:2
    line = fgetl(fid);
end
j=1;

while line ~= -1 
    %line=strrep(line,'"','');
    words = strsplit(line,',');
    humanGeneIDconversion(j,:)=[words(1),words(2)];
    words = {};
    line = {};
    line = fgetl(fid);
    j=j+1;
end
fclose(fid);