DEgenes=mouse_RNAseq_table((mouse_RNAseq_table.mouse_SCBATvsWAT_Pvalue < 0.05),:);
DEgenesHiInBAT=DEgenes((DEgenes.mouse_SCBATvsWAT_Log_Ratio > 0),:);
DEgenesMouseHuman=HumanMouse_combinedRNAseq_table((HumanMouse_combinedRNAseq_table.Mouse_SCBATvsWAT_PV < 0.05),:);
DEgenesMH2HiInBAT=DEgenesMouseHuman((DEgenesMouseHuman.Mouse_SCBATvsWAT_LR > 0),:);
plot(DEgenesMH2HiInBAT.Mouse_SCBATvsWAT_LR,DEgenesMH2HiInBAT.Human_SCBATvsWAT_LR,'b.')
hold on
xlabel('Log ratio of scBAT/WAT in mouse');
ylabel('Log ratio of scBAT/WAT in human');
