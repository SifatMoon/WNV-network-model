function Net=Net_Import3(File,N,string)


L=load(File); L1=L(:,1); L2=L(:,2); L3=L(:,3);

if strcmp(string,'directed')
[ Nei , I1 , I2 , d ] = NeighborhoodData2 ( N , L1 , L2,L3 ); I1=I1'; I2=I2';
Neigh{1}=Nei;
Net={Neigh,I1,I2,d};
end


if strcmp(string,'undirected')
dum=L1;
L1=[L1;L2];
L2=[L2;dum];
 L3=[L3;L3];
[ Nei , I1 , I2 , d ] = NeighborhoodData2 ( N , L1 , L2,L3 ); I1=I1'; I2=I2';
Neigh{1}=Nei;
Net={Neigh,I1,I2,d};
end

end