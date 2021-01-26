function [netForFlyway,l]= netForFlyway(nodes,distBwSub1,NoOfSubNet, k,flywaymap,netForFlyway,RandomNum)
% network for powerlaw-kernel biased by flyway network model
% Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus
l=1;

distBwSub=distBwSub1/10000;
distbwsumnonzerodia=distBwSub + 99999999999*eye(49,49);
xmin= min(min(distbwsumnonzerodia));
probability1 = zeros(NoOfSubNet,NoOfSubNet);
for i= 1: NoOfSubNet
    for j= i+1: NoOfSubNet
        dist= distBwSub(i,j);
        probability1(i,j) = ((k-1)/xmin)*(dist/xmin)^-k;%http://www.tandfonline.com/doi/full/10.1080/00107510500052444?scroll=top&needAccess=true
        probability1(j,i) = probability1(i,j);
    end
end

probability = probability1.*flywaymap*2;
for i =1 : size (nodes)
    for j= i : size (nodes)
        %                  prob=probability(nodes(i,2),nodes(j,2))) && (probability(nodes(i,2),nodes(j,2));
        if(nodes(i,2)~=(nodes(j,2)))
            if  (probability(nodes(i,2),nodes(j,2)) >=RandomNum(i,j))
                netForFlyway(l,1)= i;
                netForFlyway(l,2)= j;
                netForFlyway(l,3)= 0;
                
                l=l+1;
            end
            if  (probability(nodes(j,2),nodes(i,2)) >=RandomNum(j,i))
                netForFlyway(l,1)= j;
                netForFlyway(l,2)= i;
                netForFlyway(l,3)= 0;
                
                l=l+1;
                
            end
        end
    end
end
l=l-1;
