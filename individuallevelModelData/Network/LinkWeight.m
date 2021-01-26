function adjList = LinkWeight(area,time, adjList,nodes, TempByLocationCel, kk,subA2merge)
% calculate link weight of a region corresponding to the temperature of that region; links in high
% temperature region have high weights and vice versa
% Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus
s= size(adjList,1);
% listWithWeight= zeros(s, 3);
% listWithWeight(:,1:2)= adjList(:,1:2);
%     TempByLocationCel=(TempByLocation-32)*(5/9);
%     kk=.055;
minarea= 4001;
mintemp=12;
for i= 1: s
    % k1= nodes(adjList(i,1),2);
    if(nodes(adjList(i,2),3)==1)&&(nodes(adjList(i,1),3)==1)
        k1= nodes(adjList(i,1),2);
        k2= nodes(adjList(i,2),2);
        if(k1~=k2)
            temp= TempByLocationCel(k2, time);
            weight=kk/10000000000000000000;
            if(temp>mintemp)
                weight= kk*(temp-mintemp);
                %         weight=0.2;
            end
            adjList(i, 3)= weight;
        else
            temp= TempByLocationCel(k2, time);
            weight=kk/10000000000000000000;
            %                 if(k2==4)
            %                     kk=kk*0.9;
            %                 end
            if(temp>mintemp)
                weight= kk*(temp-mintemp);%*subA2merge(k2,1)/area(k2)*300;
                %         weight=0.2;
            end
            adjList(i, 3)= weight;
        end
    else
        adjList(i, 3)= 0;
    end
end