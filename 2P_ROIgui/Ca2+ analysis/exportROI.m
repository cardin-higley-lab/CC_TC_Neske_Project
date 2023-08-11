


function exportROI(ROI,outName)


[fname,pnameOut] = uiputfile('path.txt','Choose an output path...');

for i=1:length(ROI)
    ROImat(i,:)=ROI(i).fmean;
    ROIx(i)=ROI(i).centerPos(1);
    ROIy(i)=ROI(i).centerPos(2);
end
ROImat=double(ROImat');
ROIx=double(ROIx');
ROIy=double(ROIy');

dir=cd;
cd(pnameOut);
outNamex=strcat(outName,'x.txt');
outNamey=strcat(outName,'y.txt');
outNamemat=strcat(outName,'.txt');

save(outNamemat,'ROImat','-ascii');
save(outNamex,'ROIx','-ascii');
save(outNamey,'ROIy','-ascii');
cd(dir);

end