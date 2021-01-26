function llist = edgelist(adjmat)
% Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus
kk=size(adjmat);
k=kk(1,1);
llist=zeros(sum(sum(adjmat)), 2);
track=1;
for i=1:k
    for j=1:k
        if (adjmat(i,j)==1)
            llist(track,1)=j;
            llist(track,2)=i;
            track=track+1;
            
        end
    end
end
% %fprintf('%d  ', llist)
% file= fopen(filename,'w');
% %fprintf(file,'source target type\n');
% fprintf(file,'%d %d\n', llist');
% fclose(file);