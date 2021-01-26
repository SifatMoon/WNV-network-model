function k=rnd_draw(p)

% Random catagorial selection with distribution p
% Faryad Darabi Sahneh
% Kansas State University
% Last Modified: Sep 2013
% Copyright (c) 2013, Faryad Darabi Sahneh. All rights reserved. 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted

    if size(p,1)~=1
        p=p';
    end;
    a=[0,cumsum(p(1:end-1))]/sum(p);
    b=cumsum(p)/sum(p);
    toss=rand;
    k=find(a<toss & b>=toss);
end