function [NetForExpKer,l]= netForExpKer(nodes,distBwSub1,NoOfSubNet,k,NetForExpKer,RandomNum )
% network for exponential kernel network model
% Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus
l=1;

distBwSub=distBwSub1/10000;

probability = zeros(NoOfSubNet,NoOfSubNet);
for i= 1: NoOfSubNet
    for j= i+1: NoOfSubNet
        dist= distBwSub(i,j);
        probability(i,j) = k*exp(- k*dist);
        probability(j,i) = probability(i,j);
    end
end
probability = probability*2;
for i =1 : size (nodes)
    for j= i : size (nodes)
        %                  prob=probability(nodes(i,2),nodes(j,2))) && (probability(nodes(i,2),nodes(j,2));
        if(nodes(i,2)~=(nodes(j,2)))
            if  (probability(nodes(i,2),nodes(j,2)) >=RandomNum(i,j))
                NetForExpKer(l,1)= i;
                NetForExpKer(l,2)= j;
                NetForExpKer(l,3)= 0;
                
                l=l+1;
            end
            if  (probability(nodes(j,2),nodes(i,2)) >=RandomNum(j,i))
                NetForExpKer(l,1)= j;
                NetForExpKer(l,2)= i;
                NetForExpKer(l,3)= 0;
                
                l=l+1;
                
            end
        end
    end
end
l=l-1;
