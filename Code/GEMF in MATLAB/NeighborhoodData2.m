function [ Nei , I1 , I2 , d ] = NeighborhoodData2 ( N , L1 , L2,L3 )

% Find Neighborhood represention of graph from adjacency list
% Faryad Darabi Sahneh
% Kansas State University
% Last Modified: Sep 2013
% Copyright (c) 2013, Faryad Darabi Sahneh. All rights reserved. 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted
[Nodes,pr]=sort ( L1 );
Nei=[L2(pr),L3(pr)]';


l = length ( Nodes ) ;
d = zeros ( N , 1 ) ;
I1 = zeros ( N , 1 ) ;
I2 = zeros ( N , 1 ) ;
i = 1 ;
while i < l
    node = Nodes( i ) ;
    I1 ( node ) = i ;
    while Nodes( i + 1 ) == Nodes ( i )
       % node= node
        d ( node ) = d ( node ) + 1 ;
        i = i + 1 ;
        if i == l ;
            break ;
        end
    end
    i = i + 1 ;
end
if i == l
    node = Nodes( i ) ;
    I1 ( node ) = i ;
    d ( node ) = 0 ;
end
I2 = I1 + d ;
d = I2 - I1 + ( I1 ~= 0 ) ;
    
    
    