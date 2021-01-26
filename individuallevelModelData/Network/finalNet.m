function final= finalNet(NodeWithActivity, totalNet)
% Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus


l=0;
for i= 1: size(NodeWithActivity,1)
    if(NodeWithActivity(i,2)==2)
        [row, col]= find(totalNet==i) ;
        totalNet(row,:)=0;
        %         final(:,i)=0;
        %         final(i,:)=0;
        l=l+size(row,1);
    end
end
%uu= unique(totalNet(:,1));
%ll= histc(totalNet(:,1),uu);

final= zeros(size(totalNet,1)-l,2);
k=1;
for j=1: size(totalNet,1)
    if(totalNet(j,1)~=0)
        final(k,:)= totalNet(j,:);
        k=k+1;
    end
    
end