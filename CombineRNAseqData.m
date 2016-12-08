%read in Chen lab RNA-seq data from files and save in tables
%sorted VN.vs.IB_DE and SW.vs.VN_DE in excel by 1) Gene ID and 2) scBAT counts
mouseWATvsScBAT=readCSV_ChenData('SW.vs.VN_DE_excelsorted.csv');
mouseScBATvsIBAT=readCSV_ChenData('VN.vs.IB_DE.excelsorted.csv');
%running the following yields 0, which means that the two files contain the same genes
%and scBAT data in the same order and can be combined in a single table
sum(~(str2double(mouseScBATvsIBAT(2:end,4))==str2double(mouseWATvsScBAT(2:end,3))))

Entrez_ID=mouseScBATvsIBAT(2:end,1);
Gene=mouseScBATvsIBAT(2:end,2);
mouse_iBAT_mean_count=str2double(mouseScBATvsIBAT(2:end,3));
mouse_scBAT_mean_count=str2double(mouseScBATvsIBAT(2:end,4));
mouse_SCBATvsIBAT_FoldChange=str2double(mouseScBATvsIBAT(2:end,5));
%fold change is the ratio of the higher value to the lower value in this
%dataset, set to negative value if iBAT is larger; the following converts
%these values to scBAT/iBAT ratio
for i=1:length(mouse_SCBATvsIBAT_FoldChange)
    if mouse_SCBATvsIBAT_FoldChange(i) < 0
        mouse_SCBATvsIBAT_FoldChange1 = 0 - mouse_SCBATvsIBAT_FoldChange(i);
        mouse_SCBATvsIBAT_FoldChange(i) = 1./mouse_SCBATvsIBAT_FoldChange1;
    end
end
mouse_SCBATvsIBAT_Log_Ratio=log2(mouse_SCBATvsIBAT_FoldChange);
mouse_SCBATvsIBAT_Pvalue=str2double(mouseScBATvsIBAT(2:end,6));
mouse_SCBATvsIBAT_table=table(Entrez_ID,Gene,mouse_scBAT_mean_count,mouse_iBAT_mean_count,mouse_SCBATvsIBAT_Log_Ratio,mouse_SCBATvsIBAT_Pvalue);

%Entrez_ID=mouseWATvsScBAT(2:end,1);
%Gene=mouseWATvsScBAT(2:end,2);
%mouse_scBAT_mean_count=str2double(mouseWATvsScBAT(2:end,3));
mouse_WAT_mean_count=str2double(mouseWATvsScBAT(2:end,4));
mouse_WATvsSCBAT_FoldChange=str2double(mouseWATvsScBAT(2:end,5));
%fold change is the ratio of the higher value to the lower value in this
%dataset, set to negative value if scBBAT is larger; the following converts
%these values to WAT/scBAT ratio
for j=1:length(mouse_WATvsSCBAT_FoldChange)
    if mouse_WATvsSCBAT_FoldChange(j) < 0
        mouse_WATvsSCBAT_FoldChange1 = 0 - mouse_WATvsSCBAT_FoldChange(j);
        mouse_WATvsSCBAT_FoldChange(j) = 1./mouse_WATvsSCBAT_FoldChange1;
    end
end
mouse_WATvsSCBAT_Log_Ratio=log2(mouse_WATvsSCBAT_FoldChange);
%invert to get ratio of scBAT/WAT
mouse_SCBATvsWAT_Log_Ratio=0-mouse_WATvsSCBAT_Log_Ratio;
mouse_SCBATvsWAT_Pvalue=str2double(mouseWATvsScBAT(2:end,6));
mouse_SCBATvsWAT_table=table(Entrez_ID,Gene,mouse_scBAT_mean_count,mouse_WAT_mean_count,mouse_SCBATvsWAT_Log_Ratio,mouse_SCBATvsWAT_Pvalue);

%put it all together in one table
mouse_RNAseq_table=table(Entrez_ID,Gene,mouse_scBAT_mean_count,mouse_WAT_mean_count,mouse_iBAT_mean_count,mouse_SCBATvsWAT_Log_Ratio,mouse_SCBATvsWAT_Pvalue,mouse_SCBATvsIBAT_Log_Ratio,mouse_SCBATvsIBAT_Pvalue);

%read in RNA-seq data from Shinoda et al. into cell matrix humanScBATvsWAT
readShinoda;
%convert to table
RefSeq_ID=humanScBATvsWAT(:,1);
Gene={};
Gene=humanScBATvsWAT(:,2);
human_scBAT_mean_FPKM=str2double(humanScBATvsWAT(:,3));
human_WAT_mean_FPKM=str2double(humanScBATvsWAT(:,4));
human_SCBATvsWAT_Log_Ratio=str2double(humanScBATvsWAT(:,5));
human_SCBATvsWAT_table=table(RefSeq_ID,Gene,human_scBAT_mean_FPKM,human_WAT_mean_FPKM,human_SCBATvsWAT_Log_Ratio);

sortedbycount_scbat=sortrows(mouse_SCBATvsWAT_table, 3, 'descend');
%the following creates a table of the top 500 genes in scBAT by mean read count number
top500scBATonly=sortedbycount_scbat(1:500,:);
k=1;
l=1;
%the following creates a table of the top 500 genes in scBAT that are also
%expressed 2fold higher in scBAT than WAT (analogous to Shinoda data >>
%they only provided data for genes showing 2fold higher expression in
%differentiated brown adipocytes over white adipocytes
top500scBAT_2foldhigherthanWAT={};
while k<=500
    if sortedbycount_scbat{l,5}>1
        top500scBAT_2foldhigherthanWAT=[top500scBAT_2foldhigherthanWAT;sortedbycount_scbat(l,:)];
        k=k+1;
    end
    l=l+1;
end
%find human homologs of top 500 genes
top_500_Entrez_ID=str2double(top500scBATonly.Entrez_ID);
top_500_Entrez_ID_2fold=str2double(top500scBAT_2foldhigherthanWAT.Entrez_ID);
%When searching for a webtool to convert Entrez_ID to genbank accession
%numbers for the blast search, I found a tool that would not only do that,
%but would also convert between species; instead of blasting the mouse
%genbank sequences against the human genome, I opted to take advantage of
%this tool.  Part of my reasoning for doing so is because several genbank
%accession numbers were associated with each Entrez ID, and I'm not sure
%how to choose the best one.  While there are also multiple human ortholog 
%genbank ID's associated with each Entrez ID, it makes more sense to simply
%determine which of these orthologs are in the Shinoda dataset than to blast 
%all of the mouse genbank accession numbers.

%to obtain the human orthologs in the 500 gene sets, I uploaded them to the
%following website at the National Cancer Institute:
%https://biodbnet-abcc.ncifcrf.gov/db/dbOrtho.php
%using the following options: input=Gene ID; output=Genbank nucleotide
%accession; input organism: mus musculus; output organism: homo sapiens
%following is script to read in the files generated by the website:
readhumanOrth4top500file
%this script generates a cell matrix with Entrez ID in 1st cell and all of
%the Genbank accession numbers in the second cell; generate file containing
%mouse Entrez_ID, gene name, and median-normalized cpm and corresponding human Genbank
%accession number, gene name, and median-normalized FPKM

%mouseSCBATmed=median(mouse_SCBATvsIBAT_table.mouse_scBAT_mean_count,'omitnan');
%humanSCBATmed=median(human_SCBATvsWAT_table.human_scBAT_mean_FPKM,'omitnan');

%write file for top 500 scBAT only (highest counts in mouse scBAT)
fid=fopen('CountsOfTop500mouseSCBATgenesInHumanSCBAT.csv','w');
PrintString='mouse_Entrez_ID,mouse_Gene,mouse_scBAT_mean_count,mouse_WAT_mean_count,mouse_scBATvsWAT_LR,human_Genbank_accession,human_Gene,human_scBAT_mean_FPKM,human_WAT_mean_FPKM,human_scBATvsWAT_LR';
fprintf(fid,PrintString);
NumMatches=0;
NoMatches=0;

for m=1:length(top500_humanGenbank)
    %make sure everything is in same order
    inTopTable=str2double(top500scBATonly.Entrez_ID(m));
    inOrthoFile=str2double(top500_humanGenbank(m,1));
    if inTopTable==inOrthoFile
        found=0;
        words={};
        line={};
        line=strrep(top500_humanGenbank(m,2),' ','');
        line=char(line);
        line=strrep(line,';;',';');
        words=strsplit(line,';');
        %NormMouseCounts=top500scBATonly.mouse_scBAT_mean_count(m)/mouseSCBATmed;
        for n=1:length(words)
            ID2check={};
            ID2check=char(words(n));
            for o=1:length(human_SCBATvsWAT_table.RefSeq_ID)
                ScanHumanDataset=char(human_SCBATvsWAT_table.RefSeq_ID(o));
                if strcmp(ID2check,ScanHumanDataset)
                    printString={};
                    %NormHumanCounts=human_SCBATvsWAT_table.human_scBAT_mean_FPKM(o)/humanSCBATmed;
                    printString=[ '\n' char(top500scBATonly.Entrez_ID(m)) ',' char(top500scBATonly.Gene(m)) ',' num2str(top500scBATonly.mouse_scBAT_mean_count(m)) ',' num2str(top500scBATonly.mouse_WAT_mean_count(m)) ',' num2str(top500scBATonly.mouse_SCBATvsWAT_Log_Ratio(m)) ',' char(human_SCBATvsWAT_table.RefSeq_ID(o)) ',' char(human_SCBATvsWAT_table.Gene(o)) ',' num2str(human_SCBATvsWAT_table.human_scBAT_mean_FPKM(o)) ',' num2str(human_SCBATvsWAT_table.human_WAT_mean_FPKM(o)) ',' num2str(human_SCBATvsWAT_table.human_SCBATvsWAT_Log_Ratio(o))];
                    fprintf(fid,printString);
                    found=1;
                    NumMatches=NumMatches+1;
                end
            end
        end
        if found==0
            printString=[ '\n' char(top500scBATonly.Entrez_ID(m)) ',' char(top500scBATonly.Gene(m)) ',' num2str(top500scBATonly.mouse_scBAT_mean_count(m)) ',' num2str(top500scBATonly.mouse_WAT_mean_count(m)) ',' num2str(top500scBATonly.mouse_SCBATvsWAT_Log_Ratio(m)) ',human ortholog not found or not included in human RNA-seq dataset'];
            NoMatches=NoMatches+1;
            fprintf(fid,printString);
        end
    end
end
fclose(fid);

%write file for top 500 genes by count in mouse scBAT that are also 2-fold
%enriched in scBAT over WAT >> because only 2-fold enriched genes provided
%by Shinoda
fid=fopen('CountsOfTop500mouseSCBAT_2foldInHumanSCBAT.csv','w');
PrintString='mouse_Entrez_ID,mouse_Gene,mouse_scBAT_mean_count,mouse_WAT_mean_count,mouse_scBATvsWAT_LR,human_Genbank_accession,human_Gene,human_scBAT_mean_FPKM,human_WAT_mean_FPKM,human_scBATvsWAT_LR';
fprintf(fid,PrintString);
NumMatches=0;
NoMatches=0;

for m=1:length(top500_humanGenbank_2fold)
    %make sure everything is in same order
    inTopTable=str2double(top500scBAT_2foldhigherthanWAT.Entrez_ID(m));
    inOrthoFile=str2double(top500_humanGenbank_2fold(m,1));
    if inTopTable==inOrthoFile
        found=0;
        words={};
        line={};
        line=strrep(top500_humanGenbank_2fold(m,2),' ','');
        line=char(line);
        line=strrep(line,';;',';');
        words=strsplit(line,';');
        %NormMouseCounts=top500scBAT_2foldhigherthanWAT.mouse_scBAT_mean_count(m)/mouseSCBATmed;
        for n=1:length(words)
            ID2check={};
            ID2check=char(words(n));
            for o=1:length(human_SCBATvsWAT_table.RefSeq_ID)
                ScanHumanDataset=char(human_SCBATvsWAT_table.RefSeq_ID(o));
                if strcmp(ID2check,ScanHumanDataset)
                    printString={};
                    %NormHumanCounts=human_SCBATvsWAT_table.human_scBAT_mean_FPKM(o)/humanSCBATmed;
                    printString=[ '\n' char(top500scBAT_2foldhigherthanWAT.Entrez_ID(m)) ',' char(top500scBAT_2foldhigherthanWAT.Gene(m)) ',' num2str(top500scBAT_2foldhigherthanWAT.mouse_scBAT_mean_count(m)) ',' num2str(top500scBAT_2foldhigherthanWAT.mouse_WAT_mean_count(m)) ',' num2str(top500scBAT_2foldhigherthanWAT.mouse_SCBATvsWAT_Log_Ratio(m)) ',' char(human_SCBATvsWAT_table.RefSeq_ID(o)) ',' char(human_SCBATvsWAT_table.Gene(o)) ',' num2str(human_SCBATvsWAT_table.human_scBAT_mean_FPKM(o)) ',' num2str(human_SCBATvsWAT_table.human_WAT_mean_FPKM(o)) ',' num2str(human_SCBATvsWAT_table.human_SCBATvsWAT_Log_Ratio(o))];
                    fprintf(fid,printString);
                    found=1;
                    NumMatches=NumMatches+1;
                end
            end
        end
        if found==0
            printString=[ '\n' char(top500scBAT_2foldhigherthanWAT.Entrez_ID(m)) ',' char(top500scBAT_2foldhigherthanWAT.Gene(m)) ',' num2str(top500scBAT_2foldhigherthanWAT.mouse_scBAT_mean_count(m)) ',' num2str(top500scBAT_2foldhigherthanWAT.mouse_WAT_mean_count(m)) ',' num2str(top500scBAT_2foldhigherthanWAT.mouse_SCBATvsWAT_Log_Ratio(m)) ',human ortholog not found or not included in human RNA-seq dataset'];
            NoMatches=NoMatches+1;
            fprintf(fid,printString);
        end
    end
end
fclose(fid);

%the following creates a combined table of all the data for the genes
%reported by Shinoda et al. except those where it is unclear which gene in 
%mouse dataset is the ortholog or those where no ortholog was found; also writes
%file called CountsOfHumanSCBATgenesInMouseSCBAT.csv with this information

%find mouse homologs of human genes

%to obtain the mouse orthologs of genes in this set, I uploaded them to the
%following website at the National Cancer Institute:
%https://biodbnet-abcc.ncifcrf.gov/db/dbOrtho.php
%using the following options: input=Gene ID; output=Genbank nucleotide
%accession; input organism: mus musculus; output organism: homo sapiens
%following is script to read in the file generated by the website:
readmouseOrth4allHumanGenes;
%this script generates a cell matrix with human Entrez ID in 1st column and the
%mouse GeneID (if found) in the 2nd column
fid=fopen('CountsOfHumanSCBATgenesInMouseSCBAT.csv','w');
PrintString='mouse_Entrez_ID,mouse_Gene,mouse_scBAT_mean_count,mouse_WAT_mean_count,mouse_iBAT_mean_count,mouse_scBATvsWAT_Log_Ratio,mouse_scBATvsWAT_Pvalue,mouse_scBATvsiBAT_Log_Ratio,mouse_scBATvsiBAT_Pvalue,human_Genbank_accession,human_Gene,human_scBAT_mean_FPKM,human_WAT_mean_FPKM,human_scBATvsWAT_Log_Ratio';
fprintf(fid,PrintString);
NumMatches=0;
NoMatches=0;

%humanWATmed=median(human_SCBATvsWAT_table.human_WAT_mean_FPKM,'omitnan');
%mouseWATmed=median(mouse_SCBATvsWAT_table.mouse_WAT_mean_count,'omitnan');
mouseIBATmed=median(mouse_SCBATvsIBAT_table.mouse_iBAT_mean_count,'omitnan');
p=1;
for m=1:length(humanGeneIDconversion) 
    %make sure everything is in same order
    inTopTable=char(human_SCBATvsWAT_table.RefSeq_ID(m));
    inOrthoFile=char(humanGeneIDconversion(m,1));
    if strcmp(inTopTable,inOrthoFile)
        found=0;
        words={};
        line={};
        line=strrep(humanGeneIDconversion(m,2),' ','');
        line=char(line);
        line=strrep(line,';;',';');
        words=strsplit(line,';');
        %NormHumanCounts2=human_SCBATvsWAT_table.human_scBAT_mean_FPKM(m)/humanSCBATmed;
        %NormHumanWATcounts=human_SCBATvsWAT_table.human_WAT_mean_FPKM(m)/humanWATmed;
        if length(words)==1
            index=[];
            index = find(str2double(mouse_RNAseq_table.Entrez_ID)==str2double(words));
            for b=1:length(index)
                MouseEntrezID(p)=str2double(words);
                MouseGene(p)=mouse_RNAseq_table.Gene(index(b));
                MouseSCBATcount(p)=mouse_RNAseq_table.mouse_scBAT_mean_count(index(b));
                MouseWATcount(p)=mouse_RNAseq_table.mouse_WAT_mean_count(index(b));
                MouseIBATcount(p)=mouse_RNAseq_table.mouse_iBAT_mean_count(index(b));
                MouseSCBATvsWAT_LR(p)=mouse_RNAseq_table.mouse_SCBATvsWAT_Log_Ratio(index(b));
                MouseSCBATvsWAT_PV(p)=mouse_RNAseq_table.mouse_SCBATvsWAT_Pvalue(index(b));
                MouseSCBATvsIBAT_LR(p)=mouse_RNAseq_table.mouse_SCBATvsIBAT_Log_Ratio(index(b));
                MouseSCBATvsIBAT_PV(p)=mouse_RNAseq_table.mouse_SCBATvsIBAT_Pvalue(index(b));
                HumanRefSeqID(p)=human_SCBATvsWAT_table.RefSeq_ID(m);
                HumanGene(p)=human_SCBATvsWAT_table.Gene(m);
                HumanSCBATcount(p)=human_SCBATvsWAT_table.human_scBAT_mean_FPKM(m);
                HumanWATcount(p)=human_SCBATvsWAT_table.human_WAT_mean_FPKM(m);
                HumanSCBATvsWAT_LR(p)=human_SCBATvsWAT_table.human_SCBATvsWAT_Log_Ratio(m);
                printString=[ '\n' num2str(MouseEntrezID(p)) ',' char(MouseGene(p)) ',' num2str(MouseSCBATcount(p)) ',' num2str(MouseWATcount(p)) ',' num2str(MouseIBATcount(p)) ',' num2str( MouseSCBATvsWAT_LR(p)) ',' num2str(MouseSCBATvsWAT_PV(p)) ',' num2str( MouseSCBATvsIBAT_LR(p)) ',' num2str(MouseSCBATvsIBAT_PV(p)) ',' char(HumanRefSeqID(p)) ',' char(HumanGene(p)) ',' num2str(HumanSCBATcount(p)) ',' num2str(HumanWATcount(p)) ',' num2str(HumanSCBATvsWAT_LR(p)) ];
                fprintf(fid,printString);
                p=p+1;
            end
            %if more than 1 mouse gene found for RefSeq ID, check to see
            %how many are in mouse dataset; only add to table/file if 1
            %gene maps (otherwise, don't know which gene it is)
        else if length(words)>1
                howmany=(sum(str2double(mouse_RNAseq_table.Entrez_ID)==str2double(words(1)))==1);
                if howmany == 1
                    isit=words(1);
                end
                for h=2:length(words)
                    howmany= howmany + sum(str2double(mouse_RNAseq_table.Entrez_ID)==str2double(words(h))==1);
                    if howmany ==1
                        isit=words(h);
                    end
                end
                if howmany==1
                    index=[];
                    index = find(str2double(mouse_RNAseq_table.Entrez_ID)==str2double(isit));
                for c=1:length(index)
                    MouseEntrezID(p)=str2double(isit);
                    MouseGene(p)=mouse_RNAseq_table.Entrez_ID(index(c));
                    MouseSCBATcount(p)=mouse_RNAseq_table.mouse_scBAT_mean_count(index(c));
                    MouseWATcount(p)=mouse_RNAseq_table.mouse_WAT_mean_count(index(c));
                    MouseIBATcount(p)=mouse_RNAseq_table.mouse_iBAT_mean_count(index(c));
                    MouseSCBATvsWAT_LR(p)=mouse_RNAseq_table.mouse_SCBATvsWAT_Log_Ratio(index(c));
                    MouseSCBATvsWAT_PV(p)=mouse_RNAseq_table.mouse_SCBATvsWAT_Pvalue(index(c));
                    MouseSCBATvsIBAT_LR(p)=mouse_RNAseq_table.mouse_SCBATvsIBAT_Log_Ratio(index(c));
                    MouseSCBATvsIBAT_PV(p)=mouse_RNAseq_table.mouse_SCBATvsIBAT_Pvalue(index(c));
                    HumanRefSeqID(p)=human_SCBATvsWAT_table.RefSeq_ID(m);
                    HumanGene(p)=human_SCBATvsWAT_table.Gene(m);
                    HumanSCBATcount(p)=human_SCBATvsWAT_table.human_scBAT_mean_FPKM(m);
                    HumanWATcount(p)=human_SCBATvsWAT_table.human_WAT_mean_FPKM(m);
                    HumanSCBATvsWAT_LR(p)=human_SCBATvsWAT_table.human_SCBATvsWAT_Log_Ratio(m);
                    printString=[ '\n' num2str(MouseEntrezID(p)) ',' char(MouseGene(p)) ',' num2str(MouseSCBATcount(p)) ',' num2str(MouseWATcount(p)) ',' num2str(MouseIBATcount(p)) ',' num2str( MouseSCBATvsWAT_LR(p)) ',' num2str(MouseSCBATvsWAT_PV(p)) ',' num2str( MouseSCBATvsIBAT_LR(p)) ',' num2str(MouseSCBATvsIBAT_PV(p)) ',' char(HumanRefSeqID(p)) ',' char(HumanGene(p)) ',' num2str(HumanSCBATcount(p)) ',' num2str(HumanWATcount(p)) ',' num2str(HumanSCBATvsWAT_LR(p)) ];
                    fprintf(fid,printString);
                    p=p+1;
                end
                end
            end
        end
    end
end
fclose(fid);
Mouse_Entrez=MouseEntrezID';
Mouse_Gene=MouseGene';
Mouse_SCBATcount=MouseSCBATcount';
Mouse_WATcount=MouseWATcount';
Mouse_IBATcount=MouseIBATcount';
Mouse_SCBATvsWAT_LR=MouseSCBATvsWAT_LR';
Mouse_SCBATvsWAT_PV=MouseSCBATvsWAT_PV';
Mouse_SCBATvsIBAT_LR=MouseSCBATvsIBAT_LR';
Mouse_SCBATvsIBAT_PV=MouseSCBATvsIBAT_PV';
Human_RefSeqID=HumanRefSeqID';
Human_Gene=HumanGene';
Human_SCBATcount=HumanSCBATcount';
Human_WATcount=HumanWATcount';
Human_SCBATvsWAT_LR=HumanSCBATvsWAT_LR';
HumanMouse_combinedRNAseq_table=table(Mouse_Entrez,Mouse_Gene,Mouse_SCBATcount,Mouse_WATcount,Mouse_IBATcount,Mouse_SCBATvsWAT_LR,Mouse_SCBATvsWAT_PV,Mouse_SCBATvsIBAT_LR,Mouse_SCBATvsIBAT_PV,Human_RefSeqID,Human_Gene,Human_SCBATcount,Human_WATcount,Human_SCBATvsWAT_LR);
save('CombinedRNAseq_workspace.mat');