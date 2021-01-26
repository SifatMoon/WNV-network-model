function [NodeActivity, subinfection]= InitialCondition(incidence,SubNetInfo,NodeActivity,subnetinfection)
% Sifat A Moon
    c=[0;cumsum(SubNetInfo(:,1))];
    subinfection=subnetinfection;
%     for ir= 1: size(subA2merge,1)
%         p= (subA2merge(ir,1)-subA2merge(ir,time+1))/subA2merge(ir,1);
%         for jr= c(ir)+1: c(ir+1)
%             if(NodeActivity(jr,2)==1)
%                 if(rand()<p)
%                     NodeActivity(jr,2)=2;
%                 end
%             end
%         end
%     endindex= randsample(prevmodelpara2,1, true, prevWeightPower);
for i= 1:size(SubNetInfo,1)
    for j= 1:incidence(i,2)
        index= randsample(c(i)+1: c(i+1),1);
         NodeActivity(index,4)=3;
         subinfection(i,1)=subinfection(i,1)+1;
    end
end
%     for ih= 1: size(SubNetInfo,1)
%         if(incidence(ih,2)>0)
%             p= incidence(ih,2)/SubNetInfo(ih,1);
%             for j= c(ih)+1: c(ih+1)
%                 if (rand()<p)
%                     NodeActivity(j,4)=2;
%                     subinfection(ih,1)=subinfection(ih,1)+1;
%                 end
%             end
%         end
%         
%     end
%     
    