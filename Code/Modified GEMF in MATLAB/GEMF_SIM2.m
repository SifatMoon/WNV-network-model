function [ts,n_index,i_index,j_index]=GEMF_SIM2(Para,Net,x0,StopCond)

% Numerical stochastic simulation of GEMF
% Faryad Darabi Sahneh
% Kansas State University
% Last Modified: Sep 2013
% Copyright (c) 2013, Faryad Darabi Sahneh. All rights reserved. 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted

M=Para{1}; q=Para{2}; L=Para{3}; A_d=Para{4}; A_b=Para{5};

Neigh=Net{1}; I1=Net{2}; I2=Net{3}; N=length(I1);

bil=zeros(M,L);
for l=1:L
    bil(:,l)=sum(A_b(:,:,l),2);
end;

bi=cell(1,M);
for i=1:M
    temp=[];
    for l=1:L
        temp=[temp,squeeze(A_b(i,:,l))'];
    end;
    bi{i}=temp;
end;

di=sum(A_d,2);
%%----------------------------
%X0=zeros(M,N);
%for i=1:N
%    X0(x0(i),i)=1;
%end;
%X=X0;
%%in version 2 network state changes to a vector %%%%%%%%changes in v.2
X=zeros(1,N);
X=x0;
%-----------------
% Finding Nq (L by N) matrix
Nq=zeros(L,N);
%for n=1:N
%    for l=1:L
%        Nln=Neigh{l}(I1(l,n):I2(l,n));
 %       Nq(l,n)=sum(X(q(l),Nln));
 %   end;
%end;
%inversion 2 where direction of links and weight of links is considered, Nq
%is calculated differently.in version 2, Neigh{l}(1,c) where
%c=I1(l,n):I2(l,n), are the neigbours of node n in layer l that can be
%affected by it. and Neigh{l}(2,c) are the weight of the links.
for n=1:N
    for l=1:L
       if I1(l,n)~=0&&X(n)==q(l)
           Nln=Neigh{l}(1,I1(l,n):I2(l,n));
           wNln=Neigh{l}(2,I1(l,n):I2(l,n));
           for i=1:length(Nln)
                Nq(l,Nln(i))=Nq(l,Nln(i))+wNln(i);
           end;
       end;
    end;
end;

%----------------------

%Rin=di*ones(1,N).*X+bil*Nq.*X;
%Ri=sum(Rin,2);
%R=sum(Ri);
%in version 2 we do not need the above rates instesd we define a new rates
Rn=zeros(1,N);
for n=1:N
  Rn(n)=di(X(n))+bil(X(n),:)*Nq(:,n);
end
R=sum(Rn);
%_____________________________________
EventNum=StopCond{2};
RunTime=StopCond{2};

s=0; Tf=0;

while Tf<RunTime %number of events

    s=s+1;
% Event Occurance
ts(s)=-log(rand)/R;
%_______________________
% in version 2 the strategy of drawing is changed
%is=rnd_draw(Ri);
%ns=rnd_draw(Rin(is,:).*X(is,:)); % only those which are is
ns=rnd_draw(Rn);
is=X(ns);
% A_d(is,:)
% Nq(:,ns)'
% squeeze(A_b(is,:,:))
js=rnd_draw(A_d(is,:)'+bi{is}*Nq(:,ns));
%___________________________________
n_index(s)=ns;
j_index(s)=js;
i_index(s)=is;

% Updateing
% Update State
%____________________in vesion 2 updating the states changed
%X(is,ns)=0; X(js,ns)=1;
X(ns)=js;
% Updating neighbors state counters and Rates
%Ri=Ri-Rin(:,ns);
%Rin(:,ns)=di.*X(:,ns)+bil*Nq(:,ns).*X(:,ns);
%Ri=Ri+Rin(:,ns);
R=R-Rn(ns);
Rn(ns)=di(js)+bil(js,:)*Nq(:,ns);
R=R+Rn(ns);

%for l=find(q==js)
%    Nln=Neigh{l}(I1(l,ns):I2(l,ns));
%    for n=Nln
%        Nq(l,n)=Nq(l,n)+1;
%        Rin(:,n)=Rin(:,n)+bil(:,l).*X(:,n);
%    end;
%    Ri=Ri+bil(:,l).*sum(X(:,Nln),2);
%end;

%for l=find(q==is)
%    Nln=Neigh{l}(I1(l,ns):I2(l,ns));
%    for n=Nln
%        Nq(l,n)=Nq(l,n)-1;
%        Rin(:,n)=Rin(:,n)-bil(:,l).*X(:,n);
%    end;
%    Ri=Ri-bil(:,l).*sum(X(:,Nln),2);
%end;
         
for l=find(q==js)
if I1(l,ns)~=0
Nln=Neigh{l}(1,I1(l,ns):I2(l,ns));
wNln=Neigh{l}(2,I1(l,ns):I2(l,ns));
    for i=1:length(Nln)
       Nq(l,Nln(i))=Nq(l,Nln(i))+wNln(i);
       Rn(Nln(i))=Rn(Nln(i))+bil(X(Nln(i)),l)*wNln(i);
       R=R+bil(X(Nln(i)),l)*wNln(i);
    end;
end;
end;

for l=find(q==is)
if I1(l,ns)~=0
Nln=Neigh{l}(1,I1(l,ns):I2(l,ns));
wNln=Neigh{l}(2,I1(l,ns):I2(l,ns));
    for i=1:length(Nln)
        Nq(l,Nln(i))=max(0,Nq(l,Nln(i))-wNln(i));
        Rn(Nln(i))=max(0,Rn(Nln(i))-bil(X(Nln(i)),l)*wNln(i));
        R=R-bil(X(Nln(i)),l)*wNln(i);
    end;
end;
end;
%R=sum(Ri);

if R<1e-6
    break;
end;

Tf=Tf+ts(s);

end