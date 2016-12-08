% Note:  this script requires workspace from running CombineRNAseqData.m >>
% saved as CombinedRNAseq_workspace.mat

%read microarray data and Affy key into tables
brown3=readtable('GSM971719_mouseBAT3.csv','Delimiter',',');
brown2=readtable('GSM971718_mouseBAT2.csv','Delimiter',',');
brown1=readtable('GSM971717_mouseBAT1.csv','Delimiter',',');
beige3=readtable('GSM971698_mousebeige3.csv','Delimiter',',');
beige2=readtable('GSM971714_mousebeige2.csv','Delimiter',',');
beige1=readtable('GSM971716_mousebeige1.csv','Delimiter',',');
key=readtable('key2wuMicroarrayData.csv','Delimiter',',');
%replace Affy ID with Entrez ID
for i=1:length(key.ENTREZ_GENE_ID)
    brown3.Var1(strcmp(brown3.Var1,key.ID(i)))=cellstr(int2str(key.ENTREZ_GENE_ID(i)));
    brown2.Var1(strcmp(brown2.Var1,key.ID(i)))=cellstr(int2str(key.ENTREZ_GENE_ID(i)));
    brown1.Var1(strcmp(brown1.Var1,key.ID(i)))=cellstr(int2str(key.ENTREZ_GENE_ID(i)));
    beige3.Var1(strcmp(beige3.Var1,key.ID(i)))=cellstr(int2str(key.ENTREZ_GENE_ID(i)));
    beige2.Var1(strcmp(beige2.Var1,key.ID(i)))=cellstr(int2str(key.ENTREZ_GENE_ID(i)));
    beige1.Var1(strcmp(beige1.Var1,key.ID(i)))=cellstr(int2str(key.ENTREZ_GENE_ID(i)));
end

%write file containing combined mouse data (average all intensities for
%same gene in each microarray dataset)
fid=fopen('AllMouseData1.csv','w');
printString='Entrez_ID,Gene,scBAT_counts,iBAT_counts,WAT_counts,BAT1_int_Wu,BAT2_int_Wu,BAT3_int_Wu,beige1_int_Wu,beige2_int_Wu,beige3_int_Wu'; 
fprintf(fid,printString);

%sort each microarray table by Entrez ID and make sure that they contain
%the same genes in the same order
sortedBrown1=sortrows(brown1,1);
sortedBrown2=sortrows(brown2,1);
sortedBrown3=sortrows(brown3,1);
sortedBeige3=sortrows(beige3,1);
sortedBeige2=sortrows(beige2,1);
sortedBeige1=sortrows(beige1,1);

%using strcmp to compare Var1 in sortedBrown1.Var1 to Var1 in all other
%tables confirmed that they are all contain the same Entrez_IDs in the same
%order

Entrez_ID_Wu=str2double(sortedBrown1.Var1);

y=1;
for i=1:length(mouse_SCBATvsIBAT_table.Entrez_ID)
    ID=str2double(mouse_SCBATvsIBAT_table.Entrez_ID(i));
    index=[];
    index = find(Entrez_ID_Wu==ID);
    if index~=0 %only include genes in both datasets
        combinedID(y)=ID;
        combinedGene(y)=mouse_SCBATvsIBAT_table.Gene(i);
        combinedSCBATcount(y)=mouse_SCBATvsIBAT_table.mouse_scBAT_mean_count(i);
        combinedIBATcount(y)=mouse_SCBATvsIBAT_table.mouse_iBAT_mean_count(i);
        combinedWATcount(y)=mouse_SCBATvsWAT_table.mouse_WAT_mean_count(i);
        brown1int(y)=nanmean(sortedBrown1.Var2(index));
        brown2int(y)=nanmean(sortedBrown2.Var2(index));
        brown3int(y)=nanmean(sortedBrown3.Var2(index));
        beige1int(y)=nanmean(sortedBeige1.Var2(index));
        beige2int(y)=nanmean(sortedBeige2.Var2(index));
        beige3int(y)=nanmean(sortedBeige3.Var2(index));
    
        printString=[ '\n' num2str(combinedID(y)) ',' char(combinedGene(y)) ',' num2str(combinedSCBATcount(y)) ',' num2str(combinedIBATcount(y)) ',' num2str(combinedWATcount(y)) ',' num2str(brown1int(y)) ',' num2str(brown2int(y)) ',' num2str(brown3int(y)) ',' num2str(beige1int(y)) ',' num2str(beige2int(y)) ',' num2str(beige3int(y)) ];
        fprintf(fid,printString); 
        y=y+1;
        i/length(mouse_SCBATvsIBAT_table.Entrez_ID)
    end
end
fid(close);

%generate table of combined mouse gene expression data
Entrez=combinedID';
Gene=combinedGene';
scBATcounts=combinedSCBATcount';
iBATcounts=combinedIBATcount';
WATcounts=combinedWATcount';
Brown1intensity=brown1int';
Brown2intensity=brown2int';
Brown3intensity=brown3int';
Beige1intensity=beige1int';
Beige2intensity=beige2int';
Beige3intensity=beige3int';
CombinedMouseData_table=table(Entrez,Gene,scBATcounts,iBATcounts,WATcounts,Brown1intensity,Brown2intensity,Brown3intensity,Beige1intensity,Beige2intensity,Beige3intensity); 
