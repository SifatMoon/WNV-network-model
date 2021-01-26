function [gemf, T, StateCount,subinfection]= GEMF(adjList,N, x0,nodes,subinfection,dlt,lmd)
%sifat
%Susceptible-exposed-infected-recovered network model
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus
subinfection=subinfection*0;
initialInfection= size(find(x0==3));

[ Nei , I1 , I2 , d ] = NeighborhoodData2 ( N , adjList(:,1) , adjList(:,2),adjList(:,3) ); I1=I1'; I2=I2';
Neigh{1}=Nei;

Net={Neigh,I1,I2,d};


delta=6/dlt;% curing rate; dlt day to cure so rate is(1/dlt); 6 day = 1 time unit so rate will be (6/dlt)  %https://www.researchgate.net/profile/James_Hyman/publication/289125573_bergsman2016MATHEMATICAL_MODEL_FOR_THE_SPREAD_OF_WNV/links/56894d2908aebccc4e17064f/bergsman2016MATHEMATICAL-MODEL-FOR-THE-SPREAD-OF-WNV.pdf
beta=1;
lambda=6/lmd;

Para=Para_SEIR(delta,lambda,beta); M=Para{1};

runtime=5;% 1 time unit = 6 day
StopCond={'RunTime',runtime};
if(initialInfection~=0)
    [ts,n_index,i_index,j_index]=GEMF_SIM2(Para,Net,x0,StopCond);
    gemf= x0;
    timePeriod=1;
    T=[0,cumsum(ts)];
    %%% ----------post processing--------------%%%%
    for i= 1:size(n_index,2)
        %subinfection(:,i+1)=0;
        gemf(n_index(i))= j_index(i);
        if(T(i)<=timePeriod)&&(T(i)<=runtime)
            if(j_index(i)==3)
                
                subinfection(nodes(n_index(i),2), timePeriod)= subinfection(nodes(n_index(i),2),timePeriod)+1;
                
            end
        else
            while(T(i)>timePeriod)
                timePeriod= timePeriod+1;
                subinfection(:, timePeriod)= 0;
            end
            if(j_index(i)==3)
                subinfection(nodes(n_index(i),2), timePeriod)= subinfection(nodes(n_index(i),2),timePeriod)+1;
            end
        end
        
        
    end
    for tp= timePeriod+1:runtime
        subinfection(:,tp)=0;
    end
    
    [T, StateCount]=Post_Population(x0,M,N,ts,i_index,j_index);
    
else
    subinfection(:,:)=0;
    gemf=x0;
    T=30;
    StateCount(1,1)= size(find(x0==1),2);
    StateCount(2,1)= size(find(x0==2),2);
    StateCount(3,1)= size(find(x0==3),2);
    StateCount(4,1)= size(find(x0==4),2);
end



