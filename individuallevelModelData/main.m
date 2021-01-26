function main= main()
% This is the starting module
% Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus

clear all
clc

Tb1= readtable('LocationStates.xlsx'); % location of 29 states
A1=table2array(Tb1(:,2:3));
Tb2= readtable('StatePopulationARobin2016.xlsx');% American robin populaton data
A2=table2array(Tb2(:,4:15));
area= table2array(Tb2(:,2));
humanP= table2array(Tb2(:,3)); % human population data
incidence=zeros(size(humanP,1),2);

Tb3= readtable('TempState2016.xlsx'); % Average mothly temperature data
TempByLocation=table2array(Tb3(:,2:13));
Tb4= readtable('Flyway.xlsx'); % Flyway data
flyway= table2array(Tb4(:,:));
Tb5= readtable('IncidenceByState2016.xlsx'); %observed WNV human incidence data
data=table2array(Tb5(:,23+2:44+2));

TempByLocationCel=(TempByLocation-32)*(5/9); % temperature convertion: from Farenheit to celcius
time= 6; % simulation duration 6 month
[nodes,SubNetInfo,distBwSub,NoOfSubNet] =NodesInfo(A1, A2,TempByLocationCel); % preparing the nodes in the network;
flywaymap= FlyWAyMap(flyway,distBwSub,A1);
NumofTotalNode = sum(SubNetInfo(:,1));
NetAdjlist= zeros(200000,3);
%%

%incidence in south dakota and texas because incidence in 21 and 22
%%%%%%%%%%%%--------------- initial Condition--------------%%%%%%%%%%%
incidence(2,2)=3; % initial condition
incidence(40,2)=3;%initial condition

subnetinfection = zeros(size(SubNetInfo, 1),1);
subA2merge= [SubNetInfo(:,1), round(A2*0.015)+5];
subA2merge(:,1+1:time-1+1)=0;
subA2merge(:,11+1:12+1)=0;
%%
%   localNet= table2array(Tb6(:,:));
%     Tb6= readtable('localNet2016.xlsx');
%     Tb7= readtable('InitialNodeActivity2016.xlsx');
%     InitialNodeActivity=table2array(Tb7(:,:));
%     Tb8= readtable('Initialsubinfection2016.xlsx');
%     Initialsubinfection=table2array(Tb8(:,:));
%%%%%%%%%%%%------ sub-network generation -------------%%%%%%%%%%%%%%%%%
localNet = localNetER(nodes, SubNetInfo);
[InitialNodeActivity,Initialsubinfection]=InitialCondition(incidence,SubNetInfo,nodes,subnetinfection);
%%


%incidence,SubNetInfo,NodeActivity,subnetinfection
abcSMC(time, area, data, TempByLocationCel,flywaymap, NumofTotalNode, nodes, subA2merge, localNet, InitialNodeActivity, Initialsubinfection,SubNetInfo, distBwSub, NoOfSubNet, NetAdjlist)
main=0;

