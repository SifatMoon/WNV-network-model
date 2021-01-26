function Para=Para_SEIR(delta,lambda,beta)


M=4; q=[3]; L=length(q);

A_d=zeros(M); A_d(3,4)=delta; A_d(2,3)=lambda;
A_b=zeros(M,M,L); A_b(1,2,1)=beta; 


Para={M,q,L,A_d,A_b};