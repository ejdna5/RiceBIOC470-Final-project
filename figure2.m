plot ((HumanMouse_combinedRNAseq_table.Mouse_SCBATcount/median(HumanMouse_combinedRNAseq_table.Mouse_SCBATcount)),(HumanMouse_combinedRNAseq_table.Human_SCBATcount/nanmedian(HumanMouse_combinedRNAseq_table.Human_SCBATcount)),'b.')
hold on
ylim([0 40])
xlim([0 40])
xlabel('Median-normalized Mouse scBAT Count')
ylabel('Median-normalized Human scBAT Count')
hold off;
figure
plot(HumanMouse_combinedRNAseq_table.Mouse_SCBATvsWAT_LR,HumanMouse_combinedRNAseq_table.Human_SCBATvsWAT_LR,'b.');
hold on;
xlabel('Mouse scBAT/WAT log ratio')
ylabel('Human scBAT/WAT log ratio')
mousefilter=mouse_RNAseq_table(~isnan(mouse_RNAseq_table.mouse_SCBATvsWAT_Log_Ratio),:);
mousefilter2=mouse_RNAseq_table(~isinf(mousefilter.mouse_SCBATvsWAT_Log_Ratio),:);
mouse2fold=mousefilter2(mousefilter2.mouse_SCBATvsWAT_Log_Ratio>1,:);
