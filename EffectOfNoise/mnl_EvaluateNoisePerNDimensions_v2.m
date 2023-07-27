function [Dimension]=mnl_EvaluateNoisePerNDimensions_v2(maxDim)

%Function to estimate the noise generated per cell.
%Fixed Variables
NumPoints=10000;
Spreads=[0.1,0.5,1,2,4,8];
NoiseValues=0:0.01:1; %I'm goingto add the noise up to 25%

szSp=size(Spreads);
szN=size(NoiseValues,2);
%% Generate Poission Distributions For Each Set of Dimensions
for i=1:maxDim
    [Dimension(i).Cells]=mnl_GeneratePossionsNchannelsWithNoise(NumPoints,Spreads,i,NoiseValues);
    legnames{i}=sprintf('%d%s',i,' dimensions');
end
%% %For each trial calculate the distance to the 'true' value, -> Want to see histograms for example cells
[Dimension]=mnl_MeasureDistancesToCorrectValue(Dimension,szSp,szN,maxDim,NumPoints);
%% For each noise calculate the radius that covers 95% of cells
n_percent=95;
[Dimension]=mnl_CalculateRadiusForNPercent(Dimension,szSp(2),szN,maxDim,n_percent);
[Dimension]=mnl_CalculatePercentDiscriminable(Dimension,szSp(2),szN,NoiseValues,maxDim,NumPoints);
%% %Then compare the averages across the different dimensions
mnl_PlotDistancesToCorrectValuePerNoise(Dimension,Spreads,szSp,NoiseValues,szN,maxDim,NumPoints)
%% Plot some example cells with noise
mnl_Plot3Dexample(Dimension,[1,2,4],[0,0.2,0.4,0.6,0.8,1])
end
%% Subfunctions
function [Cells]=mnl_GeneratePossionsNchannelsWithNoise(NumPoints,Spreads,Ndim,NoiseValues)
% Code for creating  N dimensional poisson distributions at
% user-specified spreads, with noise added in at the NormMean level.
%
%Inputs
% NumPoints - number of cells
% Spreads - [pdf variances] e.g. [0 0.1 0.2 0.5 1 2 4 8]
% Ndim - Number of dimensions
% NoiseValues - The noise to add/subtract at the NormMean level [1 x n
% matrix]

%Outputs
% Cells=Structure containing cells and the spread information
% Structure Details - Copy Number = The specified poisson spread
%                   - XFPvals = n*Ndim 
%                   - NormXFPvals = normalised values (to max)
%                   - VecNormXFPvals= vector normalised values
%% Assign values
NumNoiseValues=size(NoiseValues,2);
rng(42) %Set the random seed
fprintf('Generating cells with %d dimensional colour space \n',Ndim);

%% Generate Poisson Matrices
maxSp=max(Spreads);
if maxSp*2>50
    MaxXFP=maxSp*2;
else
    MaxXFP=50;
end
x=linspace(0,MaxXFP,MaxXFP+1); %Theoretical maximum number of "fluorescent proteins" in the distribution
sz=size(Spreads);
for i2=1:sz(2)
    for i=1:MaxXFP
        MatlabPoisson(i,i2)=poisspdf(x(i),Spreads(i2));
    end
    legnames{i2}=sprintf('%s%d','Copy Number ',Spreads(i2));
end

%% Allocate Values for each cell
for i=1:sz(2) % For each Vector Conc
    a(:,i)=MatlabPoisson(:,i).*NumPoints;
    for j=1:Ndim %for each colour
        RealNum=0;
        Counter=1;
        for k=1:MaxXFP %Number of Copies of XFP in the cell - Theoretical Max needs to be very unlikely
            if round(MatlabPoisson(k,i)*NumPoints)>=1 %And if the value will be sufficient
                if RealNum<NumPoints %If there are more cells to label
                    NumCopies=round(MatlabPoisson(k,i)*NumPoints); %The number of times this value will appear
                    t2=ones(NumCopies,1)*k-1; %Say there are k-1 XFPs being produced, first one is zero
                    Temp(Counter:(Counter+NumCopies-1))=t2;
                    Counter=Counter+NumCopies;
                    RealNum=RealNum+NumCopies;
                end
            else
                chk=1;
            end
        end
        if Counter<=NumPoints %If at the end there is still some cells missing....
            Temp(Counter:NumPoints)=0;%then give them a zero
        elseif RealNum>NumPoints %Or if the number of cells has gone over the limit
            Temp=Temp(1:NumPoints);
        end
        Cells(i).CopyNumber=Spreads(i);
        Cells(i).XFPvals(:,j)=Temp(randperm(length(Temp))); %Shuffle the expression        
    end
    clear Temp
    Data=Cells(i).XFPvals;
    %Now normalise to max
    mxVals=max(Data,[],1);
    nData=Data./mxVals;
    Cells(i).NormXFPVals=nData;
    vnData=mnl_NormaliseVectors(nData);
    vnData(isnan(vnData))=0;
    %Now add noise per channel per cell (for 100 trials)
    nTrials=1;
    for j=1:nTrials
        for k=1:NumNoiseValues
            NoiseVal=NoiseValues(k);
            NoiseMatrix=(randn(NumPoints,Ndim))*NoiseVal; %This is the noise
            tData=vnData+NoiseMatrix;
            %move negatives to 0
            tidx=tData<0;
            tData(tidx)=0;
            %Vector Normalised
            vntData=mnl_NormaliseVectors(tData);
            %Switch NaNs to Zeros
            vntData(isnan(vntData))=0;
            Cells(i).Noise(k).Trial(j).VecNormXFPVals=vntData;
            clear data
        end
    end
    for j=1:NumNoiseValues
        Cells(i).Noise(j).NoiseAdded=NoiseValues(j);
    end
    %Store the baselines too
    Cells(i).VecNormXFPvals=vnData;
    clear vnData    
end
end
function [Dimension]=mnl_MeasureDistancesToCorrectValue(Dimension,szSp,szN,maxDim,NumPoints)
for i=1:maxDim %For each dimension
    for j=1:szSp(2) %For each spread
        TrueVals=Dimension(i).Cells(j).VecNormXFPvals; %Get the true values
        for k=1:szN %For each noise value
            nTrials=size(Dimension(i).Cells(j).Noise(k).Trial,2);
            NoiseDistances=nan(NumPoints,nTrials);
            for m=1:nTrials
                tNoiseVals=Dimension(i).Cells(j).Noise(k).Trial(m).VecNormXFPVals;
                %Measure the distance to the noise values
                [Distances]=mnl_MeasureEuclideanDistances(TrueVals,tNoiseVals);
                NoiseDistances(:,m)=Distances;
            end
            Dimension(i).Cells(j).Noise(k).NoiseValuesRaw=NoiseDistances;
            %Caluclate and save the mean, median, iqr, and sd (Global)
            Dimension(i).Cells(j).Noise(k).Global_Mean_Noise=mean(NoiseDistances(:),'omitnan');
            Dimension(i).Cells(j).Noise(k).Global_Median_Noise=median(NoiseDistances(:),'omitnan');
            Dimension(i).Cells(j).Noise(k).Global_StDev_Noise=std(NoiseDistances(:),'omitnan');
            Dimension(i).Cells(j).Noise(k).Global_IQR_Noise=iqr(NoiseDistances(:));
            %Caluclate and save the mean, median, iqr, and sd per cell
            Dimension(i).Cells(j).Noise(k).Mean_Noise=mean(NoiseDistances,2,'omitnan');
            Dimension(i).Cells(j).Noise(k).Median_Noise=median(NoiseDistances,2,'omitnan');
            Dimension(i).Cells(j).Noise(k).StDev_Noise=std(NoiseDistances,[],2,'omitnan');
            Dimension(i).Cells(j).Noise(k).IQR_Noise=iqr(NoiseDistances(:));
        end
    end
end
end
function [Distances]=mnl_MeasureEuclideanDistances(TrueVals,tNoiseVals)
%Calcluates the distances between points in a row wise manner
%Pdist2
tDistMat=pdist2(TrueVals,tNoiseVals);
Distances=diag(tDistMat);
end
function [Dimension]=mnl_CalculateRadiusForNPercent(Dimension,szSp,szN,maxDim,n_percent)
for i=1:maxDim
    for j=1:szSp
        for k=1:szN
            %Find radius containing n percent
            Distances=Dimension(i).Cells(j).Noise(k).Mean_Noise;
            Radius=prctile(Distances,n_percent);
            Dimension(i).Cells(j).Noise(k).Radius=Radius;
        end
    end
end
%Plot the relationship between noise and radius
figure("Name",'Noise to 95% Radius')
if maxDim>=8
    cmap=colormap(jet(maxDim));
else
    cmap=[0.5 0 1;0 0 1;0 1 1;0 1 0;1 1 0;1 0.5 0;1 0 0];
end
for i=1:szSp
    titleN=sprintf('%d copies per cell',num2str(round(Dimension(1).Cells(i).CopyNumber,1)));
    subplot(1,szSp,i)
    for j=1:maxDim
        legName{j}=sprintf('%d Dimensions',j);
        Xvals=nan(1,szN);
        Yvals=nan(1,szN); %These will be the radii
        for k=1:szN
            Xvals(k)=Dimension(j).Cells(i).Noise(k).NoiseAdded;
            Yvals(k)=Dimension(j).Cells(i).Noise(k).Radius;
        end
        plot(Xvals,Yvals,'color',cmap(j,:),'LineWidth',2)
        hold on
    end
    title(titleN)
    xlabel('Noise Added')
    ylabel('Radius to include upto the 95th Percentile')
    if i==szSp
        legend(legName)
    end
end

end
function mnl_PlotDistancesToCorrectValuePerNoise(Dimension,Spreads,szSp,NoiseValues,szN,maxDim,NumPoints)
% %Figure 1
% % For each noise subplot, plot the distances as boxplots for each dimension
% %for each spread
% for i=1:szSp(2)
%     fn=sprintf('Copy number = %d - Boxplots',round(Spreads(1),1));
%     figure('Name',fn)
%     for j=1:szN %Per Noise
%         tn=sprintf('Noise Level - %d',round(NoiseValues(j)));
%         NoiseBoxplotMatrix=nan(NumPoints,maxDim);
%         for k=1:maxDim
%             NoiseBoxplotMatrix(:,j)=Dimension(k).Cells(i).Noise(j).Mean_Noise;
%             ColName{k}=sprintf('%d dimensions',k);
%         end
%         subplot(round(sqrt(szN)),ceil(sqrt(szN)),j)
%         h=mnl_boxplot2(NoiseBoxplotMatrix,ColName,'Euclidean Distance','y','y');
%         title(tn)
%     end
% end
%Figure2 mean +/- sd for each noise (x)
% for each spread
%for each spread
if maxDim<=7
    cmap=[0.5 0 1;0 0 1;0 1 1;0 1 0;1 1 0;1 0.5 0;1 0 0];%bespoke 7 colour map
else
    cmap=colormap(jet(maxDim));
end

figure('Name','Effect of dimensions on noise')
for i=1:szSp(2) 
    subplot(round(sqrt(szSp(2))),ceil(sqrt(szSp(2))),i)
    for k=1:maxDim
        NoiseMatrix=nan(NumPoints,szN);
        LegName{k}=sprintf('%d dimensions',k);
        for j=1:szN %Per Noise
            NoiseMatrix(:,j)=Dimension(k).Cells(i).Noise(j).Mean_Noise;
        end
        %Mean of Means
        MeanVals=mean(NoiseMatrix,'omitnan');
        StdVals=std(NoiseMatrix,'omitnan');
        patch([NoiseValues,fliplr(NoiseValues)],[MeanVals+StdVals, fliplr(MeanVals-StdVals)],cmap(k,:),'FaceAlpha',0.5,'EdgeColor','none')
        hold on
        h(k)=plot(NoiseValues,MeanVals,'Color',cmap(k,:),'LineWidth',2);
    end
    tn=sprintf('Copy number = %s',num2str(round(Spreads(i),1)));
    title(tn)
    xlabel('Noise Added')
    ylabel('Euclidean Distance')
    legend(h,LegName)
end
end
function []=mnl_Plot3Dexample(Dimension,CopyNum,NoiseValues)
%This function is designed to plot a 3D example across the specified noise
%values, at the specified copy numbers
%Inputs
% Dimension - the dimension structure
% CopyNum - The copy numbers you want to plot
% NoiseValues - The noise values you want to plot
%
% Outputs
% Figures per copy number with one example cell

%% Initial Settings
nNoiseExamples=size(NoiseValues,2);
nCopyNumExamples=size(CopyNum,2);
nCop=size(Dimension(3).Cells,2);
for i=1:nCop
    CopyValues(i)=Dimension(3).Cells(i).CopyNumber;
end
% Calculate the 3D patch points
temp1=[linspace(0,1,101)', fliplr(linspace(0,1,101))', zeros(101,1)];
temp2=[fliplr(linspace(0,1,101))', zeros(101,1), linspace(0,1,101)'];
temp3=[zeros(101,1), linspace(0,1,101)',fliplr(linspace(0,1,101))'];
temp=[temp1;temp2;temp3];
SurfPoints=mnl_NormaliseVectors(temp);
%% Find the noise index
nAllNoise=size(Dimension(3).Cells(1).Noise,2);
AllNoiseVals=nan(nAllNoise,1);
for i=1:nAllNoise
    AllNoiseVals(i)=Dimension(3).Cells(1).Noise(i).NoiseAdded;
end
%% Make figures
for i=1:nCopyNumExamples
    cn=CopyNum(i);
    Ind=cn==CopyValues;% find the index point for cn
    % Get select a cell that is not zero
    nCells=size(Dimension(3).Cells(Ind).VecNormXFPvals,1);
    ChosenCell=round((rand(1,1)*nCells));
    GroundTruthPoint=Dimension(3).Cells(Ind).VecNormXFPvals(ChosenCell,:);
    chk=0;
    while chk==0
        if sum(GroundTruthPoint)<=0
            ChosenCell=round((rand(1,1)*nCells));
            GroundTruthPoint=Dimension(3).Cells(Ind).VecNormXFPvals(ChosenCell,:);
        else
            chk=1;
        end
    end

    fn=sprintf('Copy Number %s',num2str(round(cn,1)));
    figure('Name',fn)
    for j=1:nNoiseExamples
        subplot(round(sqrt(nNoiseExamples)),ceil(sqrt(nNoiseExamples)),j)
        %Plot Sphere
        fill3(SurfPoints(:,1),SurfPoints(:,2),SurfPoints(:,3),'k','FaceAlpha',0.1)
        hold on
        %Find Noise Index
        NoiseIndex=AllNoiseVals==NoiseValues(j);
        % Patch the surface 
        % Scatter plot of the trials
        nTrials=size(Dimension(3).Cells(1).Noise(1).Trial,2);
        NoisePlots=nan(nTrials,3);
        for k=1:nTrials
            NoisePlots(k,:)=Dimension(3).Cells(Ind).Noise(NoiseIndex).Trial(k).VecNormXFPVals(ChosenCell,:);
        end
        scatter3(NoisePlots(:,1),NoisePlots(:,2),NoisePlots(:,3),10,[0.8 0.8 0.8],'o','filled')
        hold on
        % Draw in the ground truth value
        plot3([0,GroundTruthPoint(1)],[0,GroundTruthPoint(2)],[0,GroundTruthPoint(3)],'k')
        scatter3(GroundTruthPoint(1,1),GroundTruthPoint(1,2),GroundTruthPoint(1,3),30,GroundTruthPoint,'filled')
        %set axis lims
        axis equal
        xlim([0 1])
        ylim([0 1])
        zlim([0 1])
        %Standard View Point
        view(-222,14)
        %Add title
        tn=sprintf('Noise Added %s',num2str(round(NoiseValues(j),1)));
        title(tn)
    end
end
end
function [Dimension]=mnl_CalculatePercentDiscriminable(Dimension,szSp,szN,NoiseValues,maxDim,NumPoints)
%Calculate the distance between each cell in the model, what percent of
%cells are within the radius
PcDiscrim_PerDimPerSpreadPerNoise_mean=nan(maxDim,szSp,szN);
PcDiscrim_PerDimPerSpreadPerNoise_StDev=nan(maxDim,szSp,szN);
for i=1:maxDim
    PcDiscrim_PerSpreadPerNoise_mean=nan(szSp,szN);
    PcDiscrim_PerSpreadPerNoise_StDev=nan(szSp,szN);
    for j=1:szSp
        ModelledCells=Dimension(i).Cells(j).VecNormXFPvals;
        AllDists=pdist2(ModelledCells,ModelledCells);
        PcDiscrim_PerNoise=nan(NumPoints,szN);
        for k=1:szN
            Thresh=Dimension(i).Cells(j).Noise(k).Radius;
            %Calculate Pc discriminable per row
            PcDiscrim=nan(NumPoints,1);
            for m=1:NumPoints
                idx=AllDists(m,:)>Thresh;
                PcDiscrim(m)=sum(idx)/(NumPoints-1);
            end
            Dimension(i).Cells(j).Noise(k).PercentDiscriminable=PcDiscrim;
            PcDiscrim_PerNoise(:,k)=PcDiscrim;
        end
        PcDiscrim_PerSpreadPerNoise_mean(j,:)=mean(PcDiscrim_PerNoise);
        PcDiscrim_PerSpreadPerNoise_StDev(j,:)=std(PcDiscrim_PerNoise);
    end
    PcDiscrim_PerDimPerSpreadPerNoise_mean(i,:,:)=PcDiscrim_PerSpreadPerNoise_mean;
    PcDiscrim_PerDimPerSpreadPerNoise_StDev(i,:,:)=PcDiscrim_PerSpreadPerNoise_StDev;
end
%Plot x=noise y=radius per dim?
figure('Name','Ratios Per Spread')
if maxDim<=7
    cmap=[0.5 0 1;0 0 1;0 1 1;0 1 0;1 1 0;1 0.5 0;1 0 0];
else
    cmap=colormap(jet(maxDim));
end

%Plot x=noise, y=pc discrim
for i=1:szSp
    subplot(1,szSp,i)
    for j=1:maxDim
        meanVal(1,:)=PcDiscrim_PerDimPerSpreadPerNoise_mean(j,i,:);
        stdVal(1,:)=PcDiscrim_PerDimPerSpreadPerNoise_StDev(j,i,:);
        patch([NoiseValues fliplr(NoiseValues)],[meanVal+stdVal fliplr(meanVal-stdVal)],cmap(j,:),'FaceAlpha',0.1,'EdgeAlpha',0)
        hold on
        p(j)=plot(NoiseValues,meanVal,'Color',cmap(j,:),'LineWidth',2);
        legname{j}=sprintf('%d Dimensions',j);
    end
    title(sprintf('Spread is %d',i))
    xlabel('Noise Value')
    ylabel('Percent Discriminable')
    ylim([0 1])
    if i==szSp
        legend(p,legname);
    end
end 
end