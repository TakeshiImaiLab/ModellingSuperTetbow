sz=size(MT_TeNT_all,2);
filePaths={};
Indicators={};
Time=nan(sz,1);
for i=1:sz
    filePaths{i,1}=MT_TeNT_all(i).BasicInfo.FileName;
    Indicators{i,1}=MT_TeNT_all(i).BasicInfo.IndicatorType;
    Time(i,1)=MT_TeNT_all(i).BasicInfo.MovieLength;
end

modelfunc= @(b,x) b(1).*(exp(b(2).*x))+b(3);
beta0= [0.0001,-1,0.001];
beta = nlinfit(x, y, modelfunc, beta0);