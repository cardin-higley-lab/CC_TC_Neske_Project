meanSW1 = ans.means_semiweighted(1);
meanSW2 = ans.means_semiweighted(2);

semSW1 = ans.sems_semiweighted(1);
semSW2 = ans.sems_semiweighted(2);

z = (meanSW2 - meanSW1)./sqrt((semSW2.^2 + semSW1.^2));

pValSW = 2*abs(1-normcdf(abs(z)));

workingDir = [];%mainDir/SemiweightedStats/GroupAvsGroupB/SemiweightedStruct/
comparison = 'groupAvsGroupB';

save(strcat(workingDir, comparison, '_pValSW.txt'),'pValSW','-ascii');