% Note: this script requires workspace from running CombineMouseData.m 

%write file containing combined mouse data with human data for all genes in human
%dataset that map to mouse (average all intensities for same gene in each microarray dataset)
fid=fopen('AllExpressionDataMouseLRPV.csv','w');
printString='Entrez_ID,Gene,mouse_scBAT_count,mouse_iBAT_counts,mouse_WAT_count,BAT1_int_Wu,BAT2_int_Wu,BAT3_int_Wu,beige1_int_Wu,beige2_int_Wu,beige3_int_Wu,Human_Genbank_Accession,Human_Gene,Human_scBAT_FPKM,Human_WAT_FPKM,MouseSCBATvsWAT_LR,MouseSCBATvsWAT_PV'; 
fprintf(fid,printString);

Entrez_ID_Wu=str2double(sortedBrown1.Var1);

y=1;
for i=1:length(HumanMouse_combinedRNAseq_table.Mouse_Entrez)
    ID=HumanMouse_combinedRNAseq_table.Mouse_Entrez(i);
    index=[];
    index = find(Entrez_ID_Wu==ID);
    if index~=0 %only include genes in both datasets
        mID(y)=ID;
        mGene(y)=HumanMouse_combinedRNAseq_table.Mouse_Gene(i);
        mSCBATcount(y)=HumanMouse_combinedRNAseq_table.Mouse_SCBATcount(i);
        mIBATcount(y)=HumanMouse_combinedRNAseq_table.Mouse_IBATcount(i);
        mWATcount(y)=HumanMouse_combinedRNAseq_table.Mouse_WATcount(i);
        mLR(y)=(HumanMouse_combinedRNAseq_table.Mouse_SCBATvsWAT_LR(i)>0);
        mPV(y)=HumanMouse_combinedRNAseq_table.Mouse_SCBATvsWAT_PV(i);
        hID(y)=HumanMouse_combinedRNAseq_table.Human_RefSeqID(i);
        hGene(y)=HumanMouse_combinedRNAseq_table.Human_Gene(i);
        hSCBATcount(y)=HumanMouse_combinedRNAseq_table.Human_SCBATcount(i);
        hWATcount(y)=Human_WATcount(i);
        mbrown1int(y)=nanmean(sortedBrown1.Var2(index));
        mbrown2int(y)=nanmean(sortedBrown2.Var2(index));
        mbrown3int(y)=nanmean(sortedBrown3.Var2(index));
        mbeige1int(y)=nanmean(sortedBeige1.Var2(index));
        mbeige2int(y)=nanmean(sortedBeige2.Var2(index));
        mbeige3int(y)=nanmean(sortedBeige3.Var2(index));
    
        printString=[ '\n' num2str(mID(y)) ',' char(mGene(y)) ',' num2str(mSCBATcount(y)) ',' num2str(mIBATcount(y)) ',' num2str(mWATcount(y)) ',' num2str(mbrown1int(y)) ',' num2str(mbrown2int(y)) ',' num2str(mbrown3int(y)) ',' num2str(mbeige1int(y)) ',' num2str(mbeige2int(y)) ',' num2str(mbeige3int(y)) ',' char(hID(y)) ',' char(hGene(y)) ',' num2str(hSCBATcount(y)) ',' num2str(hWATcount(y)) ',' num2str(mLR(y)) ',' num2str(mPV(y))];
        fprintf(fid,printString); 
        y=y+1;
    end
end
fid(close);

%generate table of combined mouse and human gene expression data
mouse_entrez=mID';
mouse_gene=mGene';
counts_mouseSCBAT=mSCBATcount';
counts_mouseIBAT=mIBATcount';
counts_mouseWAT=mWATcount';
intensity_mouseBrown1=mbrown1int';
intensity_mouseBrown2=mbrown2int';
intensity_mouseBrown3=mbrown3int';
intensity_mouseBeige1=mbeige1int';
intensity_mouseBeige2=mbeige2int';
intensity_mouseBeige3=mbeige3int';
human_refseq=hID';
human_gene=hGene';
counts_humanSCBAT=hSCBATcount';
counts_humanWAT=hWATcount';
LR=mLR';
PV=mPV';
CombinedMouseDataLRPV_table=table(mouse_entrez,mouse_gene,counts_mouseSCBAT,counts_mouseIBAT,counts_mouseWAT,intensity_mouseBrown1,intensity_mouseBrown2,intensity_mouseBrown3,intensity_mouseBeige1,intensity_mouseBeige2,intensity_mouseBeige3,human_refseq,human_gene,counts_humanSCBAT,counts_humanWAT,LR,PV); 
FilteredCombined=CombinedMouseDataLRPV_table((CombinedMouseDataLRPV_table.PV < 0.05),:);
DEgenes313=CombinedMouseDataLRPV_table((FilteredCombined.LR > 0),:);
save('CombinedAllExpressionLRPV_cumulativeWorkspace.mat');