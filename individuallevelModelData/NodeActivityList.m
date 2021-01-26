function NodeActivity= NodeActivityList(NodeActivity, time,subA2merge, SubNetInfo,randgemf)
%Change the activity of nodes in each month correponding to the
%population
% Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus
c=[0;cumsum(SubNetInfo(:,1))];
NodeActivity(:,3)=2;
for i= 1: size(subA2merge,1)
    %         if(subA2merge(i,time+1)==subA2merge(i,time))&&(time==4)
    %             continue;
    %         end
    %         if(subA2merge(i,time+1)<subA2merge(i,time))%  deactivate
    p= subA2merge(i,time+1)/subA2merge(i,1);
    for j= c(i)+1: c(i+1)
        %                 if(NodeActivity(j,2)==1)
        if(randgemf(j)<p)
            NodeActivity(j,3)=1;
            
        end
        %                 end
    end
    %activate
    %             q=(subA2merge(i,time+1)-subA2merge(i,time))/(subA2merge(i,1)-subA2merge(i,time));
    %             for k=  c(i)+1: c(i+1)
    %                 if(NodeActivity(k,2)==2)
    %                     if(rand()<q)
    %                         NodeActivity(k,2)=1;
    % %                         if(NodeActivity(k,3)==2)
    % %                             NodeActivity(k,3)=2;
    % %                         else
    % %                             NodeActivity(k,3)=1;
    % %                         end
    %                          NodeActivity(k,3)=1;
    %                     end
    %                 end
    %             end
    %         end
end

%          for kk=1: size(subA2merge,1)
%             check(kk)=size(find(NodeActivity(c(kk)+1: c(kk+1),2)==1),1);
%         end
%         time=time
%         check
%
%