function Para=Para_SIR(delta,beta)


M=3; q=[2]; L=length(q);

A_d=zeros(M); A_d(2,3)=delta;
A_b=zeros(M,M,L); A_b(1,2,1)=beta; 


Para={M,q,L,A_d,A_b};