function Para=Para_SIS(delta,beta)

M=2; q=[2]; L=length(q);
A_d=zeros(M); A_d(2,1)=delta;
A_b=zeros(M,M,L); A_b(1,2,1)=beta;

Para={M,q,L,A_d,A_b};

end