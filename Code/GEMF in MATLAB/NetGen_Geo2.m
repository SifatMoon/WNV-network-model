function Net=NetGen_Geo2(N,r)

% Generates random geometric graph network
% Faryad Darabi Sahneh
% Kansas State University
% Last Modified: Sep 2013
% Copyright (c) 2013, Faryad Darabi Sahneh. All rights reserved. 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted

x=rand(1,N); y=rand(1,N);

r2=r^2;

l=0;L1=[]; L2=[];
for i=1:N
    for j=i+1:N
        d2=(x(i)-x(j))^2+(y(i)-y(j))^2;
        if d2<=r2
            l=l+1;
            L1(l)=i; L2(l)=j;
            l=l+1;
            L1(l)=j; L2(l)=i;
        end;
    end;
end;

L1=L1';L2=L2';L3=ones(length(L1),1);

[ Nei , I1 , I2 , d ] = NeighborhoodData2( N , L1 , L2,L3 ); I1=I1'; I2=I2';
Neigh{1}=Nei;

Net={Neigh,I1,I2,d};

end