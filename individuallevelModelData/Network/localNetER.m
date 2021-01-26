function localNet = localNetER(nodes,SubNetInfo )
% local sub network development module
% erdos renyi random network
% Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus
%[nodes,SubNetInfo,distBwSub,NoOfSubNet] =NodesInfo();
% Nodesinfo, SubNetInfo, distBwSub,NoOfSubNet
%sizeLN= SubNetInfo(:,1)'*(SubNetInfo(:,1)-1);
localNet= zeros(size(nodes,1),3);
c=[0;cumsum(SubNetInfo(:,1))];
%ss=[0;s];
l=1;
p=10*reallog(SubNetInfo(:,1))./SubNetInfo(:,1);
for k=1:size(c)-1
    for i= c(k)+1: c(k+1)
        for j= c(k)+1: c(k+1)
            if (i~=j) %&& (nodes(i, 2)== nodes(j, 2))
                r= rand();
                if(r<=p(k))
                    localNet(l,1)= i;
                    localNet(l,2) = j;
                    l= l+1;
                    localNet(l,1)= j;
                    localNet(l,2) = i;
                    l= l+1;
                end
            end
        end
    end
end


