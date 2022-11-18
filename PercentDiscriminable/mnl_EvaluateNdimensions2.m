function [Dimension,ThresholdsPerDim]=mnl_EvaluateNdimensions2(maxDim)
Spreads=[0.1 0.2 0.5 1 2 4 6 8];
%NumPoints=(2^maxDim)-1; %If we can discriminate On and Off cells
NumPoints=10000;
szSp=size(Spreads);
%% Generate Poission Distributions For Each Set of Dimensions
for i=1:maxDim
    [Dimension(i).Cells]=mnl_GeneratePossionsNchannels(NumPoints,Spreads,i);
    legnames{i}=sprintf('%d%s',i,' dimensions');
end
%% Displaying the Barcodes for each 
szDim=size(Dimension,2);
n=1;
figure('Name','Barcodes')
cmap=magma(2^8);
colormap(cmap)
for i=1:szDim
    szCpNum=size(Dimension(i).Cells,2);
    for j=1:szCpNum
        CellNum=size(Dimension(i).Cells(j).XFPvals,1);
        TpIm=zeros(CellNum,i);
        TpIm(:,1:i)=Dimension(i).Cells(j).NormXFPvals;
        subplot(szDim,szCpNum,n)
        imagesc(TpIm,[0 1])
        n=n+1;
        if i==1
            cpn=num2str(round(Dimension(i).Cells(j).CopyNumber,1));
            tn=sprintf('%s%s','Copy Number ',cpn);
            title(tn)
        end
    end
end
%% Measure Euclidean Distances
figure('Name','Summary of Euclidean Distances')
for i=1:maxDim %For Each Dimension
    NumEuD=(NumPoints)*((NumPoints-1)/2);
    CombinedEuD_all=nan(NumEuD,szSp(2));
    for j=1:szSp(2) %For each spread
        data=Dimension(i).Cells(j).NormXFPvals;
        [EuD_all]=mnl_GroupColourEuclidean_Simplified(data); %Euclidean Distances
        CombinedEuD_all(:,j)=EuD_all;
        Dimension(i).Cells(j).EuD_all=EuD_all;
        spreadnames{j}=num2str(Spreads(j));
    end
%     subplot(1,maxDim,i)
%     mnl_boxplot(CombinedEuD_all,spreadnames,'Euclidean Distance');% Stats Graph for the spread of Euclidean Distances
%     xlabel('Number of Copies')
%     title(legnames{i})
end
%% Identify the EuD Thresholds for each Dimension
EuDThresh=0:0.01:0.5;
szEuDT=size(EuDThresh,2);
n=1;
for i=1:maxDim
    Matrix=nan(szSp(2),szEuDT);
    n1=1;
    %Allocate Matrix
    for j=1:szSp(2)
        temp=Dimension(i).Cells(j).EuD_all;
        sz=length(temp(~isnan(temp))); %the number of measurements (pairs)
        for k=1:szEuDT
            index=find(temp>EuDThresh(k));
            NumIncluded=size(index,1);
            pct=(NumIncluded/sz)*100; %Percentile rounded to 3 sig
            DimMatrix(j,k)=pct; %2D matrix - 1st dim=CopyNumber,2nd dim=Thresh EuD, Values=Percent of cells discriminated
            List(n,:)=[i EuDThresh(k) pct Dimension(i).Cells(j).CopyNumber]; %Column 1 = dim, Column 2 = Threshold, Column 3 = Percentile, Column 4 = Copy Number
            tList(n1,:)=[EuDThresh(k) pct Dimension(i).Cells(j).CopyNumber]; %Same as above but without the dimension
            n=n+1;
            n1=n1+1;
            clear index
        end
        clear temp
    end
    ThresholdsPerDim(i).Dim=i;
    ThresholdsPerDim(i).Matrix=DimMatrix;
    ThresholdsPerDim(i).List=tList;
    clear DimMatrix
    clear tList
end
%% Plot figures
% figure('Name','Plotting all Thresholds')
% cmap=colormap(jet(maxDim));
% for i=1:maxDim
%     scatter3(ThresholdsPerDim(i).List(:,1),ThresholdsPerDim(i).List(:,2),ThresholdsPerDim(i).List(:,3),'.','MarkerFaceColor',cmap(i,:))
%     hold on
% end
% xlabel('Number of Dimensions')
% ylabel('Euclidean Threshold')
% zlabel('%age of cells discriminable')
% % surface plot
% figure('Name','Surfaces')
% colormap(jet)
% for i=1:maxDim
%     [X,Y]=meshgrid(Spreads,EuDThresh);
%     Z=ThresholdsPerDim(i).Matrix';
%     C=ones(size(Z))*i;
%     s=surf(X,Y,Z,C,'FaceAlpha',0.5);
%     s.EdgeColor=cmap(i,:);
%     hold on
%     clear X Y Z C
% end
% xlabel('Copy Numbers')
% ylabel('Euclidean Threshold')
% zlabel('%age of cells discriminable')
% legend(legnames)
% figure
% for i=1:maxDim
%     subplot(2,ceil(maxDim/2),i)
%     [X,Y]=meshgrid(Spreads,EuDThresh);
%     Z=ThresholdsPerDim(i).Matrix';
%     C=ones(size(Z))*i;
%     s=surf(X,Y,Z,C);
%     s.EdgeColor=cmap(i,:);
%     xlabel('Copy Numbers')
%     ylabel('Euclidean Threshold')
%     zlabel('%age of cells discriminable')
%     tn=sprintf('%d%s',i,'Dimensions');
%     title(tn)
% end
%% Now Plot at Seven Chosen Spreads and Three Chosen Dimensions
% PossDim=[1 2 3 4 5 6 7 8];
% PossSpreads=[Spreads(1) Spreads(2) Spreads(3) Spreads(4) Spreads(5) Spreads(6) Spreads(7) Spreads(8) Spreads(9)];
% for i=1:size(PossDim,2)
%     cDim=PossDim(i);%Chosen Dimension
%     cMatrix=ThresholdsPerDim(cDim).Matrix; %Chosen Matrix
%     cmap=colormap(jet(size(PossSpreads,2)));
%     figname=sprintf('%d%s',cDim,' Dimensions - With Variable Copy Numbers');
%     figure('Name',figname)
%     for j=1:size(PossSpreads,2)
%         cSpread=PossSpreads(j);
%         %Find What row of the matrix to look for
%         id=find(Spreads==cSpread);
%         YVals=cMatrix(id,:);
%         %YVals=log10(YVals);
%         XVals=EuDThresh;
%         %Add to figure
%         plot(XVals,YVals,'Color',cmap(j,:))
%         hold on
%         lnames{j}=sprintf('%d%s',cSpread,' Copies');
%     end
%     legend(lnames)
%     xlabel('Euclidean Distance Threshold')
%     ylabel('Percent Discriminable')
%     clear cMatrix YVals XVals lnames
% end
%% Seperate plots of discriminable vs. EuThresh
linecolours=mnl_CreateRGBdistributions(7);
linecolours=linecolours/255;
Loop=ceil(szSp(2)/7);
tLineColours=[];
for i=1:Loop
    tLineColours=[tLineColours;linecolours];
end
linecolours=tLineColours;
clear tLineColours
figure
Xvals=EuDThresh;
for i=1:maxDim
    tn=sprintf('%d%s',i,' Dimensions');
    subplot(ceil(sqrt(maxDim)),ceil(sqrt(maxDim)),i)
    for j=1:szSp(2) %For each Spread
        Yvals=ThresholdsPerDim(i).Matrix(j,:);
        plot(Xvals,Yvals,'Color',linecolours(j,:));
        hold on
        lname{j}=sprintf('%d%s',Spreads(j),' Copies');
    end
    legend(lname,'Location','eastoutside')
    title(tn)
end
end
function [ChannelColours]=mnl_CreateRGBdistributions(numChan)
C1=[255 0 0];
C2=[0 255 0;255 0 0];
C3=[0 0 255;0 255 0;255 0 0];
C4=[125 0 255;0 0 255;0 255 0;255 0 0];
C5=[125 0 255;0 0 255;0 255 255;0 255 0;255 0 0];
C6=[125 0 255;0 0 255;0 255 255;0 255 0;255 255 0;255 0 0];
C7=[125 0 255;0 0 255;0 255 255;0 255 0;255 255 0;255 125 0;255 0 0];
if numChan==1
    ChannelColours=C1;
elseif numChan==2
    ChannelColours=C2;
elseif numChan==3
    ChannelColours=C3;
elseif numChan==4
    ChannelColours=C4;
elseif numChan==5
    ChannelColours=C5;
elseif numChan==6
    ChannelColours=C6;
elseif numChan==7
    ChannelColours=C7;
end
end