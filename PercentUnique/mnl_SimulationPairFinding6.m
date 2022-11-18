function [DimensionSummary]=mnl_SimulationPairFinding6(MaxNumPoints,nDim,NumSim)
%Spreads=[0.1 0.2 0.5 1 2 4 6 8];
%Spreads=[0.2 0.5 1 2 4 6 8];
%Spreads=[2 4 6 8];
Spreads=[0.5 1 2 4];
%% Create colour map
[cmap]=mnl_CalculateColourMap(nDim);
%Spreads=[2]; %Down to 2 Copy Numbers
szSp=size(Spreads);
%% Generate Poission Distributions For Each Set of Dimensions
%PointRange=1:1:MaxNumPoints; If you want to do all cell points
PointRange=0:5:MaxNumPoints; %Every 5
PointRange(1)=1;
[Dim]=mnl_GeneratePoissons(nDim,NumSim,PointRange,Spreads);
%% Do this with EuD thresholds
disp('Now Per EuThresh')
c=1;
EuDThresh=0:0.05:0.5;
nCellsPerTrial=PointRange;
nGroups=size(nCellsPerTrial,2);
TotalNum=nGroups*nDim*NumSim*szSp(2);
st=tic;
for d=1:nDim
    for i=1:szSp(2) %For Each Spread
        for j=1:NumSim %For Each Sim
            for k1=1:nGroups %For each maximum number of cells
                k=nCellsPerTrial(k1);
                Matrix=Dim(d).Sim(j).MaxCellNum(k).Cells(i).NormXFPvals;
                [~,PercentageUnique,~]=mnl_FindPercentUnique(Matrix,EuDThresh,k);
                %Dim(d).Sim(j).MaxCellNum(k).Cells(i).EuDmatrix=EuD;
                %Dim(d).Sim(j).MaxCellNum(k).Cells(i).Pairs=Pairs;
                Dim(d).Sim(j).MaxCellNum(k).Cells(i).PercentageUnique=PercentageUnique;
                mnl_InsertProgressTrackerInLoops(c,TotalNum)
                c=c+1;
            end
        end
    end
end

%% Make a plot with the 95% CI and mean
x=PointRange;
for j=1:size(EuDThresh,2) %For Each EuD Threshold
    tThresh=round(EuDThresh(j),2);
    tThresh=num2str(tThresh);
    fn=sprintf('%s%s','Percentage Unique at EuD Thresh ',tThresh);
    figure('Name',fn) 
    for i=1:szSp(2) % For Each Spread
        SpreadNum=Spreads(i);
        subplot(1,szSp(2),i)
        if i==szSp(2)
            legnum=1;
        end
        for d=1:nDim
            DimensionSummary(d).Spread(i).NumberOfCells=PointRange;
            DimensionSummary(d).Spread(i).NumSims=NumSim;
            DimensionSummary(d).Spread(i).CopyNum=SpreadNum;
            SimCollections=nan(NumSim,MaxNumPoints);
            for k=1:NumSim %For each simulation
                %Collate the requisit Percentage Values
                for m1=1:nGroups
                    m=PointRange(m1);
                    SimCollections(k,m)=Dim(d).Sim(k).MaxCellNum(m).Cells(i).PercentageUnique(j);
                    SimCollections2(k,m1)=Dim(d).Sim(k).MaxCellNum(m).Cells(i).PercentageUnique(j);
                end
            end
            DimensionSummary(d).Spread(i).EuDThresh(j).EuValue=EuDThresh(j);
            DimensionSummary(d).Spread(i).SimCollections(j).Matrix=SimCollections;
            DimensionSummary(d).Spread(i).Mean(j,:)=nanmean(SimCollections);
            DimensionSummary(d).Spread(i).StandardDev(j,:)=nanstd(SimCollections);
            DimensionSummary(d).Spread(i).EuThresh(j).PercentUnique=SimCollections;
            yMean=nanmean(SimCollections2);
            yStd=nanstd(SimCollections2);
            %Now plot
            P2=patch([x fliplr(x)], [yMean+yStd fliplr(yMean-yStd)], cmap(d,:),'EdgeColor','none');
            hold on
            P2.FaceAlpha=0.2;
            pId(d)=plot(x,yMean,'Color',cmap(d,:),'LineWidth',2);
            if i==szSp(2)
                legnames{legnum}=sprintf('%d%s',d,' XFPs');
                %legnames{legnum+1}=sprintf('%d%s',d,' XFPs');
                legnum=legnum+1;
            end
        end
        xlim([0 MaxNumPoints])
        ylim([0 100])
        xlabel('Number of Cells')
        ylabel('Percentage Unique Per Trial')
        SpN=round(SpreadNum,1);
        SpN=num2str(SpN);
        subpTitle=sprintf('%s%s',SpN,' Copies');
        title(subpTitle)
        if i==szSp(2)
            legend(pId,legnames)
        end
    end
%     mfn=sprintf('%s%s',fn,'.fig');
%     %efn=sprintf('%s%s',fn,'.eps');
%     h=gcf;
%     savefig(h,mfn)
%     %mnl_ExportEPSdense(h,efn)
end
%% Per EuThresh
%[cmap]=mnl_GenerateShuffledColourmap(nGroups);
[cmap]=colormap(jet(nGroups));
nThresh=size(EuDThresh,2);
clear SimCollections
for d=1:nDim
    clear legnames
    fn=sprintf('%s%d%s','EuThresh_X_',d,' XFPs');
    figure('Name',fn)
    for i=1:szSp(2) % For Each Spread
        SpreadNum=Spreads(i);
        subplot(1,szSp(2),i)
        if i==szSp(2)
            legnum=1;
        end
        for j1=1:nGroups %For each number of cells
            j=PointRange(j1);
            SimCollections=nan(NumSim,nThresh);
            for k=1:NumSim %For each simulation
                %Collate the requisit Percentage Values for each EuThresh
                for m=1:nThresh
                    SimCollections(k,m)=Dim(d).Sim(k).MaxCellNum(j).Cells(i).PercentageUnique(m);
                end
            end
            yMean=nanmean(SimCollections);
            yStd=nanstd(SimCollections);
            %Now Plot
            P2=patch([EuDThresh fliplr(EuDThresh)], [yMean+yStd fliplr(yMean-yStd)], cmap(j1,:),'EdgeColor','none');
            hold on
            P2.FaceAlpha=0.2;
            pId(j1)=plot(EuDThresh,yMean,'Color',cmap(j1,:),'LineWidth',2);
            if i==szSp(2)
                legnames{legnum}=sprintf('%d%s',j,' cells');
                legnum=legnum+1;
            end
        end
        xlim([0 EuDThresh(nThresh)])
        ylim([0 100])
        xlabel('Euclidean Distance Threshold')
        ylabel('Percentage Unique Per Trial')
        SpN=round(SpreadNum,1);
        SpN=num2str(SpN);
        subpTitle=sprintf('%s%s',SpN,' Copies');
        title(subpTitle)
        if i==szSp(2)
            legend(pId,legnames)
        end
    end
end
%% Per Number of Cells
%[cmap]=mnl_GenerateShuffledColourmap(nThresh);
[cmap]=colormap(jet(nThresh));
nThresh=size(EuDThresh,2);
for d=1:nDim
    clear legnames
    clear pId
    fn=sprintf('%s%d%s','Ncells_X_',d,' XFPs');
    figure('Name',fn)
    for i=1:szSp(2) % For Each Spread
        SpreadNum=Spreads(i);
        subplot(1,szSp(2),i)
        if i==szSp(2)
            legnum=1;
        end
        for j=1:nThresh %For each Euclidean Distance
            SimCollections=nan(NumSim,nGroups);
            for k=1:NumSim %For each simulation
                %Collate the requisit Percentage Values for each EuThresh
                for m1=1:nGroups
                    m=PointRange(m1);
                    SimCollections(k,m1)=Dim(d).Sim(k).MaxCellNum(m).Cells(i).PercentageUnique(j);
                end
            end
            yMean=nanmean(SimCollections);
            yStd=nanstd(SimCollections);
            %Now Plot
            P2=patch([PointRange fliplr(PointRange)], [yMean+yStd fliplr(yMean-yStd)], cmap(j,:),'EdgeColor','none');
            hold on
            P2.FaceAlpha=0.2;
            pId(j)=plot(PointRange,yMean,'Color',cmap(j,:),'LineWidth',2);
            if i==szSp(2)
                legnames{legnum}=sprintf('%s%s','Euclidean Distance ',num2str(round(EuDThresh(j),2)));
                legnum=legnum+1;
            end
        end
        xlim([0 PointRange(nGroups)])
        ylim([0 100])
        xlabel('Euclidean Distance Threshold')
        ylabel('Percentage Unique Per Trial')
        SpN=round(SpreadNum,1);
        SpN=num2str(SpN);
        subpTitle=sprintf('%s%s',SpN,' Copies');
        title(subpTitle)
        if i==szSp(2)
            legend(pId,legnames)
        end
    end
end
end
function [cmap]=mnl_CalculateColourMap(nDim)
if nDim<=7
    cmap=nan(nDim,3);
    for i=1:nDim
        if i==1
            cmap(i,:)=[0.5 0 1];
        elseif i==2
            cmap(i,:)=[0 0 1];
        elseif i==3
            cmap(i,:)=[0 0.5 1];
        elseif i==4
            cmap(i,:)=[0 1 0];
        elseif i==5
            cmap(i,:)=[1 1 0];
        elseif i==6
            cmap(i,:)=[1 0.5 0];
        elseif i==7
            cmap(i,:)=[1 0 0];
        end
    end
else
    cmap=colormap(jet(nDim));
end
end
function [Dim]=mnl_GeneratePoissons(nDim,NumSim,PointRange,Spreads)
disp('Generating All the Poissons')
szPR=size(PointRange,2);
totalNum=nDim*NumSim*szPR;
c=1;
Dim=struct('Sim',[]);
for d=1:nDim
    for i=1:NumSim
        for j=1:szPR
            numCells=PointRange(j);
            [Dim(d).Sim(i).MaxCellNum(numCells).Cells]=mnl_GeneratePossionsNchannels(numCells,Spreads,d);
            mnl_InsertProgressTrackerInLoops(c,totalNum)
            c=c+1;
        end
    end
end
end
function [pairs,PercentageUnique,EuD]=mnl_FindPercentUnique(Matrix,EuDThresh,k)
[EuD]=mnl_GroupEuclidean_Matrixv2(Matrix,Matrix);
%Are there values below the following thresholds?
pairs=zeros(k,size(EuDThresh,2));
for m=1:k % For each cell
    for n=1:size(EuDThresh,2)
        index=find(EuD(m,:)<=EuDThresh(n));
        szI=size(index,2);
        pairs(m,n)=pairs(m,n)+szI; %row-cell# col-EuD thresh value-number of pairs
        clear index
    end
end

%Now find out the percentage of unique ones
for n=1:size(EuDThresh,2)
    index2=find(pairs(:,n)==1); %It matches with itself so that is one match
    NumUnique(n)=size(index2,1);
end
PercentageUnique=(NumUnique/k)*100;
end

function [EuD]=mnl_MeasureEuclideanDistances(Matrix1,Matrix2)
%Function to measure the Euclidean Distance between two matrices
% Input
% Matrix1/2 - rows=each trace, cols=each channel
%
% Output
% EuD - Matrix of all the euclidean distances
sz1=size(Matrix1);
sz2=size(Matrix2);
EuD=nan(sz1(1),sz2(1));
for i=1:sz1(1)
    Pos1=Matrix1(i,:);
    for j=1:sz2(1)
        Pos2=Matrix2(j,:);
        EuVal=mnl_EuclideanDistance(Pos1,Pos2);
        EuD(i,j)=EuVal;
    end
end
end
function EuD=mnl_EuclideanDistance(Pos1,Pos2)
Diffs=Pos1-Pos2;
DiffsSq=Diffs.^2;
Tot=sum(DiffsSq);
EuD=sqrt(Tot);
end