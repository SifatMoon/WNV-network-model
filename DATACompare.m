function dataCompare= DATACompare()
Tb1= readtable('IncidenceByState2016.xlsx');
A1=table2array(Tb1(:,23:44));
Tb2=readtable('ModelResult.xlsx');
A2=table2array(Tb2(:,1:26));

% Tb3= readtable('StatePopulationARobin.xlsx');
% % A2=table2array(Tb2(:,4:15));
% humanP= table2array(Tb3(:,3));
% humanPo=zeros(size(humanP,1),2);
% humanpo(:,1)= humanP;
ita=1/15;%%%******************parameter
Modeldata=zeros(size(A2,1), size(A1,2));
%Modeldata(:,1:7)=A2(:,2:8);
timeperiod=6;
datadist=zeros(size(A2,1),size(A2,2)*timeperiod);
for i= 2:size(A2,2)
    datadist(:,(i-2)*timeperiod+1)=A2(:,i)./6;
    for j= 1: timeperiod-1
        datadist(:,(i-2)*timeperiod+1+j)=datadist(:,(i-2)*timeperiod+1);
    end
end
for i= 1: size(A1,2)
    Modeldata(:,i)=sum(datadist(:,((i-1)*7+1):(i-1)*7+7),2);
end


mse= zeros(size(A2,1),1);
n=size(A1,2);
for  i= 1: size(A2,1)
   mse(i,1)= sum(((A1(i,:)-Modeldata(i,:)).^2*ita))/n;

end
plot_distfloat(mse)

function [] = plot_distfloat(Data)

    dist = zeros(1,ceil(max(Data))+1);
    
    for i=1:length(Data)
        dist(ceil(Data(i))) = dist(ceil(Data(i))) + 1;
    end
   figure
    plot(0:ceil(max(Data)),dist,'r*')