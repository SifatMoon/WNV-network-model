function Net=NetCmbn2(NetSet)

% Combines simple graphs into a multilayer graph
% Faryad Darabi Sahneh
% Kansas State University
% Last Modified: Sep 2013
% Copyright (c) 2013, Faryad Darabi Sahneh. All rights reserved. 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted

L=length(NetSet);

for l=1:L
    Neigh{l}=NetSet{l}{1}{1};
    I1(l,:)=NetSet{l}{2};
    I2(l,:)=NetSet{l}{3};
    d(:,l)=NetSet{l}{4};
  
end;

Net={Neigh,I1,I2,d};

end
    
