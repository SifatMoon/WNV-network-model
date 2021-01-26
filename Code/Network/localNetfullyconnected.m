function localNet = localNetfullyconnected(nodes,SubNetInfo )
% Sifat A Moon
% adj list
%[nodes,SubNetInfo,distBwSub,NoOfSubNet] =NodesInfo();
% Nodesinfo, SubNetInfo, distBwSub,NoOfSubNet
sizeLN= SubNetInfo(:,1)'*(SubNetInfo(:,1)-1);
localNet= zeros(sizeLN,2);
s= cumsum(SubNetInfo(:,1));
ss=[0;s];
l=1;
for i= 1: size(nodes)
    
    localNet(l:l-2+SubNetInfo(nodes(i,2),1),1)=i;
    localNet(l:l-2+SubNetInfo(nodes(i,2),1),2)= nodes(ss(nodes(i,2)) +1:ss(nodes(i,2)+1)~=i,1);
    %localNet(l:l-2+SubNetInfo(nodes(i,2),1),3)= 0;
    l= l-1+SubNetInfo(nodes(i,2),1);
end


