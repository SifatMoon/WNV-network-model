function [simulationResult,isCandidateData,totalStateCountPerDayAlliterReturn]= simulation(area,numOfTotalNode,time,m, localNet, nodes,subNetInfo,subA2Merge,distBwSub,noOfSubNet,initialNodeActivity, initialsubinfection,tempByLocationCel,netPara,betaCons,delta, lambda,flywayMap,netAdjList)
% main simulation module
% Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus
nfl=50; % number of iteration
lastMonth=10;
simulationResult=zeros(size(subNetInfo, 1),(lastMonth-(time-1))*5+1);
result=zeros(size(subNetInfo, 1),(lastMonth-(time-1))*5+1);
%%%%%-----network development---%%%
for fl=1:nfl
    randomNum = rand(numOfTotalNode,numOfTotalNode);
    if(m==1)
        [netForKer,lnet]= netForExpKer(nodes,distBwSub,noOfSubNet, netPara,netAdjList,randomNum);
    elseif(m==2)
        [netForKer,lnet]= netForPowerLaw(nodes,distBwSub,noOfSubNet, netPara,netAdjList,randomNum);
    else
        [netForKer,lnet]= netForFlyway(nodes,distBwSub,noOfSubNet, netPara,flywayMap,netAdjList,randomNum);
    end
    isCandidateData=0;
    %%%%% -------Temporal network behavior--------%%%%%%%%%%%%
    dayBegin=(time-1)*30-1;
    totalStateCountPerDayAlliterReturn= zeros(4,lastMonth*30 -dayBegin);
    if(lnet>5000)% if the number of link is for distance kernel is more than 5000
        isCandidateData=1;
        totalNet= [localNet;netForKer(1:lnet, :)];
        iter=60;
        
        
        totalstateCountPerDayalliter=zeros(4,lastMonth*30 -dayBegin);
        totalStateCountPerDayAlliterReturn= zeros(4,lastMonth*30 -dayBegin);
        randgemf=rand(numOfTotalNode,iter*(lastMonth-time+1));
        totalsubinfectionFull=zeros(size(subNetInfo, 1),(lastMonth-(time-1))*5+1);
        for i= 1: iter
            subinfectionFull=initialsubinfection;
            subinfection=initialsubinfection;
            Tfull=(time-1)*5;
            totalstateCountPerDay=zeros(4,lastMonth*30 -dayBegin);
            
            StateCountFull=[numOfTotalNode;0;5;0];
            NodeActivity=initialNodeActivity;
            for t= time:lastMonth
                totalNet1=totalNet;
                if((t==time)&&(m==3))
                    totalNet1(:,1)=totalNet(:,2);
                    totalNet1(:,2)=totalNet(:,1);
                end
                
                NodeWithActivity = NodeActivityList(NodeActivity, t,subA2Merge, subNetInfo,randgemf(:,(i -1)*(lastMonth-time+1)+t-time+1));
                listWithWeight = LinkWeight(area,t, totalNet1,NodeWithActivity, tempByLocationCel, betaCons, subA2Merge);
                % calling GEMF module
                [gemf, T, StateCount,subinfection]=GEMF(listWithWeight, numOfTotalNode, NodeActivity(:,4)', nodes, subinfection,delta,lambda);
                %post processing
                Tfull=[Tfull, (t-1)*5+T];
                subinfectionFull=[subinfectionFull,subinfection];
                StateCountFull=[StateCountFull,StateCount];
                subinfection(:,1)=subinfection(:,size(subinfection,2));
                NodeActivity(:,4)=gemf';
            end
            % %post processing
            day=(time-1)*30+1;
            stateCountPerDay=zeros(4,(lastMonth-time)*30);
            Tfull=Tfull*6;
            stateCountPerDay(:,1)=StateCountFull(:,1);
            for k=1: size(Tfull,2)
                while(Tfull(k)>day)
                    if(day>lastMonth*30)
                        break;
                    end
                    stateCountPerDay(:,day-dayBegin)=StateCountFull(:,k);
                    day=day+1;
                end
                
            end
            if(day~=lastMonth*30+1)
                for iday=day:lastMonth*30
                    stateCountPerDay(:,iday-dayBegin)=stateCountPerDay(:,day - 1 - dayBegin) ;
                end
            end
            l=size(stateCountPerDay);
            if(l(2)==151)
                totalstateCountPerDay =totalstateCountPerDay+stateCountPerDay;
            end
            totalsubinfectionFull= totalsubinfectionFull + subinfectionFull;
            totalstateCountPerDayalliter= totalstateCountPerDay+ totalstateCountPerDayalliter;
            clear  subinfectionFull;
            clear StateCountFull;
            clear T;
            clear gemf;
            clear listWithWeight;
            clear NodeWithActivity;
        end
        result=result + totalsubinfectionFull/iter;
        totalStateCountPerDayAlliterReturn = totalstateCountPerDayalliter/iter;
        
    else
        simulationResult=0;
    end
    clear netForKer;
    clear totalNet;
    clear randomNum;
    clear randgemf;
end

simulationResult= result/nfl;
clear  result;
clear totalsubinfectionFull;


