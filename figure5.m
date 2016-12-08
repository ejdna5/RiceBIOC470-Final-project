 X=table2array(CombinedAllData_table(:,[3:11,14,15]));
 Xmed=nanmedian(X);
 for i=1:length(Xmed)
    Xnorm(:,i)=X(:,i)/Xmed(i);
 end
[coeff,sc,eig]=pca(Xnorm');
labels = {'mouseSCBAT','mouseIBAT','mouseWAT','brown1','brown2','brown3','beige1','beige2','beige3','humanSCBAT','humanWAT'}; 
figure
plot(sc(:,1),sc(:,2),'b.','MarkerSize',18);
hold on
text(sc(:,1),sc(:,2),labels);
xlabel('PC1');
ylabel('PC2');
figure
plot(sc(:,1),sc(:,3),'b.','MarkerSize',18);
hold on
text(sc(:,1),sc(:,3),labels);
xlabel('PC1');
ylabel('PC3');
figure
plot(sc(:,2),sc(:,3),'b.','MarkerSize',18);
hold on
text(sc(:,2),sc(:,3),labels);
xlabel('PC2');
ylabel('PC3');