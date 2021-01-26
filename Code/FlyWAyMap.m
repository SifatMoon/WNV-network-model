
function flywaymap= FlyWAyMap(flyway,distBwSub,A1)
% develop flywaymap
% Sifat Afroj Moon
%https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006875&rev=2#sec024
%A spatio-temporal individual-based network framework for West Nile virus in the USA: Spreading pattern of West Nile virus
flywaymap= zeros(length(flyway),length(flyway));
for i=1:4
    for j= 1: length(flyway)
        for k=j+1: length(flyway)
            if(flyway(j,i)*flyway(k,i)>0)
                if((distBwSub(j,k))>=500*10^3 )
                    if   A1(j)>A1(k)
                        flywaymap(j,k)= flyway(j,i);
                    else
                        flywaymap(k,j)= flyway(k,i);
                    end
                end
            end
        end
    end
end
for i=1:length(flyway)
    for j= i+1:length(flyway)
        if((distBwSub(i,j))<=500*10^3 ) % connection in the flyway between tw long distance locations/states
            flywaymap(i,j)= 1;
            flywaymap(j,i)= 1;
        end
    end
end

flywaymap(42,35)= 1;