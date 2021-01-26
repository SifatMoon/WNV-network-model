function abcSMC= abcSMC(time, area, data, tempByLocationCel,flywayMap, numOfTotalNode, nodes, subA2Merge, localNet, initialNodeActivity, initialsubinfection, subNetInfo, distBwSub, noOfSubNet, netAdjlist)
%Approximate Bayesian Computation based on sequentianl Monte Carlo sampling
%module(ABC-SMC): parameter estimation and model selection
%Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus

%%%----------------Input for ABC-SMC method------------%%
epsilon=[4000,3000,2500, 2000, 1900, 1850, 1800, 1700,1600,1500,1220,1200,1180,1150,1120,1100,1080,1050,1020,1000 ]; %tolerance array
epsilonTime= length(epsilon);
NPerEpsilon=1000; % number of particles
NoOfModel=3;
noOfParaModel1=4; % number of parameter in the first model (exponetial kernel network model)
noOfParaModel2=4; % number of parameter in the second model (power-law kernel network model)
noOfParaModel3=4; % number of parameter in the third model (powe-law kernel biased by flyway network model)

modelPara1=NPerEpsilon;
modelPara2=NPerEpsilon;
modelPara3=NPerEpsilon;

newWeightExp= zeros(NPerEpsilon,1)+1; % store weight of  the accepted particles of the first model
newWeightPower= zeros(NPerEpsilon,1)+1;  % store weight of  the accepted particles of the second model
newWeightFlyway= zeros(NPerEpsilon,1)+1; % store weight of  the accepted particles of the third model

SMCRound=2;

for population = SMCRound: epsilonTime % SMC round
    Tb9= readtable('powerKernelPopulation.xlsx'); % accepted parameters for first model for (1: SMCRound-1)
    Tb10= readtable('expKernelPopulation.xlsx');% accepted parameters for second model for (1: SMCRound-1)
    Tb11= readtable('flywaykernelPopulation.xlsx');% accepted parameters for third model for (1: SMCRound-1)
    
    paraExpAll= zeros(NPerEpsilon, epsilonTime*noOfParaModel1);
    paraPowerAll= zeros(NPerEpsilon, epsilonTime*noOfParaModel2);
    ParaFlywayAll= zeros(NPerEpsilon, epsilonTime*noOfParaModel3);
    
    powerKernelPopulation=table2array(Tb9(:,:));
    expKernelPopulation=table2array(Tb10(:,:));
    flywaykernelPopulation=table2array(Tb11(:,:));
    
    paraExpAll(1:NPerEpsilon,1:noOfParaModel1*(population-1))=expKernelPopulation(1:NPerEpsilon,1:noOfParaModel1*(population-1));
    paraPowerAll(1:NPerEpsilon,1:noOfParaModel2*(population-1))=powerKernelPopulation(1:NPerEpsilon,1:noOfParaModel2*(population-1));
    ParaFlywayAll(1:NPerEpsilon,1:noOfParaModel3*(population-1))=flywaykernelPopulation(1:NPerEpsilon,1:noOfParaModel3*(population-1));
    
    prevModelPara1=modelPara1;
    prevModelPara2=modelPara2;
    prevModelPara3=modelPara3;
    
    expPerPop=[260,412,550,150,172,57,103,230,42,0,0,0,0,0,0,0,0,0]; % Accepted particles of first model in each SMC rounds
    powerPerPop=[1000,417,560,112,113,82,133,182,28,0,0,0,0,0,0,0,0,0,0];% Accepted particles of second model in each SMC rounds
    flyPerPop=[900,445, 621,131, 200,175,172, 466,87,];% Accepted particles of third model in each SMC rounds
    
    strecordExptb= readtable('strecordExp.xlsx');
    strecordExp= table2array(strecordExptb(:,:));
    
    strecordPowertb= readtable('strecordPower.xlsx');
    strecordPower= table2array(strecordPowertb(:,:));
    
    strecordflytb= readtable('strecordFly.xlsx');
    strecordfly= table2array(strecordflytb(:,:));
    newWeightFlyway=strecordfly(1:flyPerPop(1),1);
    
    prevWeightExp= strecordExp(1:expPerPop(1),1)/ sum(strecordExp(1:expPerPop(1),1));
    prevWeightPower= strecordPower(1:powerPerPop(1),1)/sum(strecordPower(1:powerPerPop(1),1));
    prevWeightFlyway= strecordfly(1:flyPerPop(1),1)/sum(strecordfly(1:flyPerPop(1),1));
    
    
    newWeightExp= newWeightExp*0;
    newWeightPower=newWeightPower* 0;
    newWeightFlyway= newWeightFlyway*0;
    
    prevModelPara1=expPerPop(population-1);
    prevModelPara2=powerPerPop(population-1);
    prevModelPara3=flyPerPop(population-1);
    if(population>2)
        
        %%%%prevWeight, prevPara, para, noPrevPara, noPara, model
        for jj= 3:population
            model=1;
            newWeightExp=newWeightExp*0;
            parfor ii=1:expPerPop(jj-1)
                parameter=paraExpAll(ii, (jj-2)*noOfParaModel1+1:(jj-2)*noOfParaModel1+noOfParaModel1);
                newWeightExp(ii)= calculateWeight(prevWeightExp,paraExpAll(:,(jj-3)*noOfParaModel1+1:(jj-3)*noOfParaModel1+noOfParaModel1),parameter,expPerPop(jj-2),noOfParaModel1,model);
            end
            %%
            newWeightExp(1:expPerPop(jj-1),1)=newWeightExp(1:expPerPop(jj-1),1).*strecordExp(1:expPerPop(jj-1),jj-1);
            prevWeightExp(1:expPerPop(jj-1),1)= newWeightExp(1:expPerPop(jj-1),1)/sum(newWeightExp(1:expPerPop(jj-1),1));
        end
        
        for jj= 3:population
            model=2;
            newWeightPower=newWeightPower*0;
            parfor ii=1:powerPerPop(jj-1)
                parameter=paraPowerAll(ii, (jj-2)*noOfParaModel2+1:(jj-2)*noOfParaModel2+noOfParaModel2);
                newWeightPower(ii)= calculateWeight(prevWeightPower,paraPowerAll(:,(jj-3)*noOfParaModel2+1:(jj-3)*noOfParaModel2+noOfParaModel2),parameter,powerPerPop(jj-2),noOfParaModel2,model);
            end
            %%
            newWeightPower(1:powerPerPop(jj-1),1)=newWeightPower(1:powerPerPop(jj-1),1).*strecordPower(1:powerPerPop(jj-1),jj-1);
            prevWeightPower(1:powerPerPop(jj-1),1)= newWeightPower(1:powerPerPop(jj-1),1)/sum(newWeightPower(1:powerPerPop(jj-1),1));
        end
        for jj= 3:population
            model=3;
            newWeightFlyway=newWeightFlyway*0;
            parfor ii=1:flyPerPop(jj-1)
                parameter=ParaFlywayAll(ii, (jj-2)*noOfParaModel3+1:(jj-2)*noOfParaModel3+noOfParaModel3);
                newWeightFlyway(ii)= calculateWeight(prevWeightFlyway,ParaFlywayAll(:,(jj-3)*noOfParaModel3+1:(jj-3)*noOfParaModel3+noOfParaModel3),parameter,flyPerPop(jj-2),noOfParaModel3,model);
            end
            %%
            newWeightFlyway(1:flyPerPop(jj-1),1)=newWeightFlyway(1:flyPerPop(jj-1),1).*strecordfly(1:flyPerPop(jj-1),jj-1);
            prevWeightFlyway(1:flyPerPop(jj-1),1)= newWeightFlyway(1:flyPerPop(jj-1),1)/sum(newWeightFlyway(1:flyPerPop(jj-1),1));
        end
        
    end
    for i=1:10000
        sob(area,numOfTotalNode,data,epsilon,prevWeightExp,prevWeightPower,prevWeightFlyway,paraExpAll,paraPowerAll,ParaFlywayAll,population,NoOfModel,prevModelPara1,prevModelPara2, prevModelPara3, noOfParaModel1, noOfParaModel2, noOfParaModel3,time,localNet,nodes,subNetInfo,subA2Merge,distBwSub,noOfSubNet,initialNodeActivity, initialsubinfection,tempByLocationCel,flywayMap,netAdjlist);
    end
end
abcSMC=0;

function netParaExpPrior= netParaExpPrior()
netParaExpPrior= 0.1+0.4*rand();

function netParaPowerAlphaPrior= netParaPowerAlphaPrior()
netParaPowerAlphaPrior= 2+ 2*rand();

function netParaFlywayPrior=netParaFlywayPrior()
netParaFlywayPrior=2+ 2*rand();

function betaConsPrior= betaConsPrior()
betaConsPrior= 15*rand();

function deltaPrior= deltaPrior()
deltaPrior= 4+rand();

function lambdaPrior= lambdaPrior()
lambdaPrior=25*rand();

function uniKernel= uniKernel(alpha)
uniKernel= alpha *(-1 +2*rand());

%%%%%%%%%%

function densityExp= densityExp(k)
densityExp = unifpdf(k,0.1,.5);

function densityPowerAlpha= densityPowerAlpha(k)
densityPowerAlpha = unifpdf(k,2,4);

function densityFlyway= densityFlyway(k)
densityFlyway = unifpdf(k,2,4);

function densityBeta= densityBeta(k)
densityBeta= unifpdf(k,0,15);

function densityDelta=  densityDelta(k)
densityDelta= unifpdf(k,4,5);

function densitylambda=  densitylambda(k)
densitylambda= unifpdf(k,0,25);

function desityUniKernel= desityUniKernel(k,alpha)
desityUniKernel=unifpdf(k, -alpha, alpha);
%%%%%%%%%%%%%%%%%%

function [checkData,minDistVal,proCons,minpcModelData]= checkData(data, candidateData, epsilon)
% compare the simulation data and the observed incidence data
minDistVal=0;
modelData=zeros(size(candidateData,1), size(data,2));
timePeriod=6;
dataDist=zeros(size(candidateData,1),size(candidateData,2)*timePeriod);
for i= 2:size(candidateData,2)
    dataDist(:,(i-2)*timePeriod+1)=candidateData(:,i)./timePeriod;
    for j= 1: timePeriod-1
        dataDist(:,(i-2)*timePeriod+1+j)=dataDist(:,(i-2)*timePeriod+1);
    end
end
for i= 1: size(data,2)
    modelData(:,i)=sum(dataDist(:,((i-1)*7+1):(i-1)*7+7),2);
end
proCons=1;
minpcModelData=zeros(49,22);
minpcModelData=modelData;

minDistVal=sqrt(sum(sum((cumsum(data,2)- cumsum(modelData,2)).^2)));
for i=1:100
    pc=0.1+50*rand();
    pcModelData1=poissrnd(pc*modelData);
    for ii=1:9
        pcModelData1=pcModelData1+poissrnd(pc*modelData);
    end
    pcModelData=pcModelData1/10;
    
    distVal=sqrt(sum(sum((cumsum(data,2)- cumsum(pcModelData,2)).^2)));
    if(distVal<minDistVal)
        minpcModelData=pcModelData;
        minDistVal=distVal;
        proCons=pc;
    end
end

if(minDistVal<epsilon)
    checkData=1;
end

function calculateWeight=calculateWeight(prevWeight, prevPara, para, noPrevPara, noPara, model)
%CAlculate weight for the accepted particles
vals=zeros(noPrevPara,1);
for i= 1: noPrevPara
    vals(i)= prevWeight(i);
    for j= 1: noPara
        diff= para(j)- prevPara(i,j);
        if (model ==1)
            if(j==1)
                vals(i)= vals(i)* desityUniKernel(diff, 0.06);
            elseif(j==2)
                vals(i)= vals(i)* desityUniKernel(diff, 6);
            elseif(j==3)
                vals(i)= vals(i)* desityUniKernel(diff, 0.5);
            else
                vals(i)= vals(i)* desityUniKernel(diff, 5);
            end
        elseif(model==2)
            
            if (j==1)
                vals(i)= vals(i)* desityUniKernel(diff, 0.3*2);
            elseif (j==2)
                vals(i)= vals(i)* desityUniKernel(diff, 6);
            elseif (j==3)
                vals(i)= vals(i)* desityUniKernel(diff, 0.5);
            else
                vals(i)= vals(i)* desityUniKernel(diff, 5);
            end
        else
            if(j==1)
                vals(i)= vals(i)* desityUniKernel(diff,1);
            elseif(j==2)
                vals(i)= vals(i)* desityUniKernel(diff, 8);
            elseif(j==3)
                vals(i)= vals(i)* desityUniKernel(diff, 0.5);
            else
                vals(i)= vals(i)* desityUniKernel(diff, 6);
            end
        end
    end
end

priorDensityForCurrent=0;
if (model ==1)
    priorDensityForCurrent= densityExp(para(1,1))*densityBeta(para(1,2))*densityDelta(para(1,3))*densitylambda(para(1,4));
elseif(model==2)
    priorDensityForCurrent=densityPowerAlpha(para(1,1))*densityBeta(para(1,2))*densityDelta(para(1,3))*densitylambda(para(1,4));
else
    priorDensityForCurrent=densityFlyway(para(1,1))*densityBeta(para(1,2))*densityDelta(para(1,3))*densitylambda(para(1,4));
end
calculateWeight= priorDensityForCurrent/sum(vals);




function sobv= sob(area,numOfTotalNode,data,epsilon,prevWeightExp,prevWeightPower,prevWeightFlyway,paraExpAll,paraPowerAll,ParaFlywayAll,population,NoOfModel,prevModelPara1,prevModelPara2, prevModelPara3, noOfParaModel1, noOfParaModel2, noOfParaModel3,time,localNet,nodes,subNetInfo,subA2Merge,distBwSub,noOfSubNet,initialNodeActivity, initialsubinfection,tempByLocationCel,flywayMap,netAdjlist)
modelPara1=0;
modelPara2=0;
modelPara3=0;
perticle=0;
trial=0;
paratoChange=1;
iter=1;
ParaFlywayStoW= zeros(1000,1);
ParaExpwayStoW=zeros(1000,1);
ParaPowerwayStoW=zeros(1000,1);

parameter=zeros(1,4);
while perticle<10000
    trial=trial+1;
    model=randi(NoOfModel,1);
    %         model=3;
    
    if(model==1)
        if (population==1)
            
            parameter(1,1)= netParaExpPrior();
            parameter(1,2)= betaConsPrior();
            parameter(1,3)= deltaPrior();
            parameter(1,4)= lambdaPrior();
            
        else
            index= randsample(prevModelPara1,1, true, prevWeightExp(1:prevModelPara1));
            parameter(1,1:noOfParaModel1) = paraExpAll(index,((population-2)*noOfParaModel1+1):((population-2)*noOfParaModel1+noOfParaModel1));
            paratoChange= randi(noOfParaModel1,1);
            
            if (paratoChange ==1)
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(.07);%alpha
            elseif (paratoChange ==2)
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(3);%beta
            elseif (paratoChange ==3)
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(0.25);%Delta
            else
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(3);%lambda
            end
            
        end
        if (densityExp(parameter(1,1))*densityBeta(parameter(1,2))*densityDelta(parameter(1,3))*densitylambda(parameter(1,4)))>0
            icd=0;
            for cd=1:iter
                [candidateData,isCandidateData,totalStateCountPerDay]= simulation(area,numOfTotalNode,time,model,localNet,nodes,subNetInfo,subA2Merge,distBwSub,noOfSubNet,initialNodeActivity, initialsubinfection,tempByLocationCel,parameter(1,1),parameter(1,2),parameter(1,3),parameter(1,4),flywayMap,netAdjlist);
                if(isCandidateData==1)
                    [checkData1,minDistVal,proCons, modeldata]= checkData(data, candidateData, epsilon(population));
                    
                    if (checkData1==1)&&(length(find(sum(candidateData,2)>0))>=30)
                        ss=length(find(sum(candidateData,2)>0))
                        icd=icd+1;
                    end
                end
            end
            
            
            if(icd>0)
                perticle=perticle+1;
                modelPara1=modelPara1+1;
                thmodeldata((modelPara1-1)*49 +1:(modelPara1-1)*49 +49,1:22)=modeldata;
                thmodeldata1(1:49,modelPara1)=sum(modeldata,2);
                expProCons(modelPara1,1)=proCons;
                % alltotalstateCountPerDay((modelPara1-1)*4 +1:(modelPara1-1)*4 +4,1:151)= totalStateCountPerDay;
                paraExpAll(modelPara1,((population-1)* noOfParaModel1+1):((population-1)* noOfParaModel1+noOfParaModel1))=parameter(1,1:noOfParaModel1) ;
                ParaExpwayStoW(modelPara1,1)=icd/iter;
                
                
            end
        end
    elseif(model ==2)
        if (population==1)
            parameter(1,1)= netParaPowerAlphaPrior();
            parameter(1,2)= betaConsPrior();
            parameter(1,3)= deltaPrior();
            parameter(1,4)= lambdaPrior();
            
        else
            index= randsample(prevModelPara2,1, true, prevWeightPower(1:prevModelPara2));
            parameter(1,1:noOfParaModel2) = paraPowerAll(index,((population-2)*noOfParaModel2+1):((population-2)*noOfParaModel2+noOfParaModel2));
            paratoChange= randi(noOfParaModel2,1);
            
            if (paratoChange ==1)
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(1);%alpha
            elseif (paratoChange ==2)
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(3);%beta
            elseif (paratoChange ==3)
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(0.25);%Delta
            else
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(2);%lambda
            end
        end
        
        
        
        if (densityPowerAlpha(parameter(1,1))*densityBeta(parameter(1,2))*densityDelta(parameter(1,3))*densitylambda(parameter(1,4)))>0
            
            icd=0;
            for cd=1:iter
                [candidateData,isCandidateData,totalStateCountPerDay]= simulation(area,numOfTotalNode,time,model,localNet,nodes,subNetInfo,subA2Merge,distBwSub,noOfSubNet,initialNodeActivity, initialsubinfection,tempByLocationCel,parameter(1,1),parameter(1,2),parameter(1,3),parameter(1,4),flywayMap,netAdjlist);
                if(isCandidateData==1)
                    [checkData1,minDistVal,proCons, modeldata]= checkData(data, candidateData, epsilon(population));
                    
                    if (checkData1==1)&&(length(find(sum(candidateData,2)>0))>=30)
                        ss=length(find(sum(candidateData,2)>0))
                        icd=icd+1;
                    end
                end
            end
            
            
            if(icd>0)
                perticle=perticle+1;
                modelPara2=modelPara2+1;
                thmodeldata((modelPara2-1)*49 +1:(modelPara2-1)*49 +49,1:22)=modeldata;
                thmodeldata2(1:49,modelPara2)=sum(modeldata,2);
                PoproCons(modelPara2,1)=proCons;
                %alltotalstateCountPerDay((modelPara2-1)*4 +1:(modelPara2-1)*4 +4,1:151)= totalStateCountPerDay;
                paraPowerAll(modelPara2,((population-1)* noOfParaModel2+1):((population-1)* noOfParaModel2+noOfParaModel2))=parameter(1,1:noOfParaModel2) ;
                ParaPowerwayStoW(modelPara2,1)=icd/iter;
                
                
            end
        end
    else
        if (population==1)
            
            parameter(1,1)= netParaFlywayPrior();
            parameter(1,2)= betaConsPrior();
            parameter(1,3)= deltaPrior();
            parameter(1,4)= lambdaPrior();
            
        else
            index= randsample(prevModelPara3,1, true, prevWeightFlyway(1:prevModelPara3));
            parameter(1,1:noOfParaModel3) = ParaFlywayAll(index,((population-2)*noOfParaModel3+1):((population-2)*noOfParaModel3+noOfParaModel3));
            paratoChange= randi(noOfParaModel2,1);
            
            if (paratoChange ==1)
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(1);%alpha
            elseif (paratoChange ==2)
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(2);%beta
            elseif (paratoChange ==3)
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(0.25);%Delta
            else
                parameter(1,paratoChange)= parameter(1,paratoChange)+ uniKernel(2);%lambda
            end
            
        end
        
        if (densityFlyway(parameter(1,1))*densityBeta(parameter(1,2))*densityDelta(parameter(1,3))*densitylambda(parameter(1,4)))>0
            %                          parpool(2);
            
            icd=0;
            for cd=1:iter
                [candidateData,isCandidateData,totalStateCountPerDay]= simulation(area,numOfTotalNode,time,model,localNet,nodes,subNetInfo,subA2Merge,distBwSub,noOfSubNet,initialNodeActivity, initialsubinfection,tempByLocationCel,parameter(1,1),parameter(1,2),parameter(1,3),parameter(1,4),flywayMap,netAdjlist);
                if(isCandidateData==1)
                    [checkData1,minDistVal,proCons, modeldata]= checkData(data, candidateData, epsilon(population));
                    
                    if (checkData1==1)&&(length(find(sum(candidateData,2)>0))>=30)
                        ss=length(find(sum(candidateData,2)>0))
                        icd=icd+1;
                        minDistVal=minDistVal
                    end
                end
            end
            
            if(icd>0)
                perticle=perticle+1;
                modelPara3=modelPara3+1;
                thmodeldata((modelPara3-1)*49 +1:(modelPara3-1)*49 +49,1:22)=modeldata;
                thmodeldata3(1:49,modelPara3)=sum(modeldata,2);
                flyproCons(modelPara3,1)=proCons;
                ParaFlywayAll(modelPara3,((population-1)* noOfParaModel3+1):((population-1)* noOfParaModel3+noOfParaModel3))=parameter(1,1:noOfParaModel3) ;
                ParaFlywayStoW(modelPara3,1)=icd/iter;
            end
        end
    end
end
filename=  ['Exp',num2str(trial),num2str(10*rand()),'.txt'];
file= fopen(filename,'w');
%fprintf(file,'source target type\n');
expprint= [paraExpAll(:,((population-1)* noOfParaModel1+1):((population-1)* noOfParaModel1+noOfParaModel1)), ParaExpwayStoW];
fprintf(file,'%f %f %f %f %f\n', expprint');
fclose(file);

filename=  ['Power',num2str(trial),num2str(10*rand()),'.txt'];
file= fopen(filename,'w');
%fprintf(file,'source target type\n');
powerprint= [paraPowerAll(:,((population-1)* noOfParaModel2+1):((population-1)* noOfParaModel2+noOfParaModel2)), ParaPowerStoW];
fprintf(file,'%f %f %f %f %f\n', powerprint');
fclose(file);

filename=  ['Flyway',num2str(trial),num2str(10*rand()),'.txt'];
file= fopen(filename,'w');
%fprintf(file,'source target type\n');
flyprint= [ParaFlywayAll(:,((population-1)* noOfParaModel3+1):((population-1)* noOfParaModel3+noOfParaModel3)), ParaFlywayStoW];
fprintf(file,'%f %f %f %f %f\n', flyprint');
fclose(file);

sobv=0;



