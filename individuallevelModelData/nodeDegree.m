function []= nodeDegree()
% Sifat Afroj Moon
    Tb1= readtable('totalNet.xlsx');
A1=table2array(Tb1(:,1:2));
m= max(max(A1));
adjmat= zeros(m,m);
for i= 1: size(A1,1)
    adjmat(A1(i,1), A1(i,2))=1;
 
end
nodeDegree= sum(adjmat,2);
plot_dist(nodeDegree);
function [] = plot_dist(Data)

    dist = zeros(1,max(Data)+1);

    for i=1:length(Data)
        dist(Data(i)+1) = dist(Data(i)+1) + 1;
    end
   figure
    bar(0:max(Data),dist/sum(dist))