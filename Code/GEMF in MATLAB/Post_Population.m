function [T, StateCount]=Post_Population(x0,M,N,ts,i_index,j_index)

% Postprocessing of stochastic simulation. From events to population size
% Faryad Darabi Sahneh
% Kansas State University
% Last Modified: Sep 2013
% Copyright (c) 2013, Faryad Darabi Sahneh. All rights reserved. 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted

X0=zeros(M,N);
for i=1:N
    X0(x0(i),i)=1;
end;

T=[0,cumsum(ts)];
StateCount=zeros(M,length(ts)+1);
StateCount(:,1)=sum(X0,2);
for k=1:length(ts)
    DX=zeros(M,1); DX(i_index(k))=-1; DX(j_index(k))=1;
    StateCount(:,k+1)=StateCount(:,k)+DX;
end;

end
