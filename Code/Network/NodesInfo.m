function [Nodesinfo, SubNetInfo, distBwSub,NoOfSubNet]= NodesInfo(A1,A2, TempByLocationCel)
%% preparing module for nodes and subnetworks
% Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus

A2=A2*0.015 +5; % American population scaling
% A2=A2*0.0015 +5;
A2= round(A2);
NoOfSubNet= size(A1,1);
for i= 1:49
    mx=0;
    for j= 6:10
        if (TempByLocationCel(i, j)>=18)
            if(A2(i,j)>mx)
                mx=A2(i,j);
            end
        end
    end
    NodesinSub(i,1)=mx;
end

% NodesinSub=round((max(A2(:,6:8), [], 2)*4)/100)+5;% rescaled american robin population divided by 1000%10 is margin which i took so no subnode will not lost
SubNetInfo =[NodesinSub, A1];
ToNode= sum(NodesinSub);
Nodesinfo=zeros(ToNode,4);
%NodeActivity(1:ToNode,1)=1:ToNode;
Nodesinfo(:,3)=1;%1= active ; 2= deactive
Nodesinfo(:,4)=1;%1= susceptible in susceptible-exposed-infected-recovered (SEIR) epidemic model
j=1;
for i = 1: size(NodesinSub, 1)
    Nodesinfo(j:j-1+NodesinSub(i),1)=j:j-1+NodesinSub(i);%node number
    Nodesinfo(j:j-1+NodesinSub(i),2)=i; % subnet number
    j=j+ NodesinSub(i);
end
distBwSub = zeros(NoOfSubNet, NoOfSubNet);
for i= 1: NoOfSubNet
    for j= i+1:NoOfSubNet
        distBwSub(i,j) = DistBwn2GeoCordinate(A1(i, 1), A1(i, 2), A1(j, 1), A1(j,2));
        distBwSub(j,i)= distBwSub(i,j);
    end
end