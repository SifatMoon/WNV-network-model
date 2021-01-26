%function [gemf, T, StateCount,subinfection]= GEMF(adjList,N, x0,nodes,subinfection)
    %sifat
    N=379;
    adjList=load('edgesweight.txt');
    %Network generation
    x0=ones(1,N); x0(1:10)=2; %initial condition of the network, "1" represents susceptible "2" represents infected

    initialInfection= size(find(x0==2));
    
    [ Nei , I1 , I2 , d ] = NeighborhoodData2 ( N , adjList(:,1) , adjList(:,2),adjList(:,3) ); I1=I1'; I2=I2';
    Neigh{1}=Nei;
    
    Net={Neigh,I1,I2,d};
   
    delta=1/10;
    beta=.09;
    Para=Para_SIR(delta,beta); M=Para{1};
%     lambda1=10;
% delta=1; beta=3/lambda1;
% Para=Para_SIS(delta,beta); M=Para{1};
  
    
%     tic;
    StopCond={'RunTime',100};
    if(initialInfection~=0)
        [ts,n_index,i_index,j_index]=GEMF_SIM2(Para,Net,x0,StopCond);
        
%         toc;
        gemf= x0;
      
%         for i= 1:size(n_index,2)
%             %subinfection(:,i+1)=0;
%             gemf(n_index(i))= j_index(i);
%             if(j_index(i)==2)
%                subinfection(nodes(n_index(i),2))= subinfection(nodes(n_index(i),2))+1;
%             else
%                 subinfection(nodes(n_index(i),2))= subinfection(nodes(n_index(i),2))-1;
%             end
%             
%         end
%         Post Processing
      StatesPlot=[1,2,3];
        [T, StateCount]=Post_Population(x0,M,N,ts,i_index,j_index);
plot(T,StateCount(StatesPlot,:)/N);
    else
        gemf=x0;
        T=30;
        StateCount(1,1)= size(find(x0==1),2);
        StateCount(2,1)= size(find(x0==2),2);
        StateCount(3,1)= size(find(x0==3),2);
    end
    
    
    
    