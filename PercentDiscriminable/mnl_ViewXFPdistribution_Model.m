function []=mnl_ViewXFPdistribution_Model(Dimension)
NumDim=size(Dimension,2);
figure
for i=1:NumDim
    NumCopies=size(Dimension(i).Cells,2);
    for j=1:NumCopies
        VecNormMatrix=Dimension(i).Cells(j).NormXFPvals;
        tn=sprintf('%d%s%d',i,' Dimensions, Copy Number - ',Dimension(i).Cells(j).CopyNumber);
        subplot(NumDim,NumCopies,((NumCopies*(i-1))+j))
        mnl_CumulativePlotMatrix(VecNormMatrix)
        title(tn)
        xlabel('Vector Normalised Value')
    end
end
%% Merge the dimensions
figure('Name','Spread of the XFP values')
for i=1:NumCopies
    CopyNum=Dimension(1).Cells(i).CopyNumber;
    CopyNum=round(CopyNum,1);
    CopyNum=num2str(CopyNum);
    CopyTitle=sprintf('%s%s',CopyNum,' Copy Numbers');
    for j=1:NumDim
        XFPvals=Dimension(j).Cells(i).XFPvals;
        %Plot the XFP values
        subplot(2,NumCopies,i)
        mnl_CumulativePlotMatrix2(XFPvals)
        title(CopyTitle)
        ylim([0 100])
        %Plot the Norm Vec Values
        VecNormVals=Dimension(j).Cells(i).NormXFPvals;
        subplot(2,NumCopies,i+NumCopies)
        mnl_CumulativePlotMatrix2(VecNormVals)
        title(CopyTitle)
        ylim([0 100])
    end
end
%% Do it per dimension
figure('Name','Per Dimension')
for i=1:NumDim
    subplot(1,NumDim,i)
    DimTitle=sprintf('%d%s',i,' Dimensions');
    for j=1:NumCopies
        CopyNum=Dimension(i).Cells(j).CopyNumber;
        CopyNum=round(CopyNum,1);
        CopyNum=num2str(CopyNum);
        CopyTitle=sprintf('%s%s',CopyNum,' Copy Numbers');
        legnames{j}=CopyTitle;
        %Plot the Norm Vec Values
        VecNormVals=Dimension(i).Cells(j).NormXFPvals;
        mnl_CumulativePlotMatrix2(VecNormVals)
        hold on
    end
    title(DimTitle)
    ylim([0 100])
    xlim([0 1])    
    if i==NumDim
        legend(legnames)
    end
end
%% Barcodes
figure
for i=1:NumDim
    for j=1:NumCopies
        RGB=zeros(50,i,3);
        [Map]=mnl_SelectRGBcolours(i);
        VecNormMatrix=Dimension(i).Cells(j).NormXFPvals;
        for k=1:50
            for m=1:i
                RGB(k,m,:)=VecNormMatrix(k,m).*Map(m,:);
            end
        end
        subplot(NumDim,NumCopies,(NumCopies*(i-1))+j)
        image(RGB)
        tn=sprintf('%s%d','Copy Numbers - ',Dimension(i).Cells(j).CopyNumber)
        title(tn)
        clear RGB
    end
end
        
%%
%3D
prompt='Do you want 3D Plots? - y/n';
ThreeD_Plots=input(prompt,'s');
if strcmp('y',ThreeD_Plots)==1
    i=5;
    Cmap=Dimension(3).Cells(i).NormXFPvals;
    figure
    for j=1:10000
        scatter3(Dimension(3).Cells(i).NormXFPvals(j,1),Dimension(3).Cells(i).NormXFPvals(j,2),Dimension(3).Cells(i).NormXFPvals(j,3),100,'.','MarkerFaceColor',Cmap(j,:),'MarkerEdgeColor',Cmap(j,:))
        hold on
    end
    tn=sprintf('%s%d','Copy Number',Dimension(3).Cells(i).CopyNumber);
    title(tn)
    view(140,40)
end
%Ternary
prompt='Do you want Ternary Plots? - y/n';
TP_Plots=input(prompt,'s');
if strcmp(TP_Plots,'y')==1
    for i=1:NumCopies
        figure
        [h,hg,htick]=terplot;
        x=Dimension(3).Cells(i).NormXFPvals(:,1);
        y=Dimension(3).Cells(i).NormXFPvals(:,2);
        z=Dimension(3).Cells(i).NormXFPvals(:,3);
        [Group]=Convert2Ratios([x y z]);
        hter=ternaryc(Group(:,1),Group(:,2),Group(:,3));
        for j=1:10000
            set(hter(j),'marker','o','markeredgecolor','none','markerfacecolor',[x(j) y(j) z(j)],'markersize',6)
        end
        % Now create labels
        terlabel('XFP 1','XFP 2','XFP 3');
        tn=sprintf('%s%d','Copy Number',Dimension(3).Cells(i).CopyNumber);
        title(tn)
    end
end
end
%% Sub function
function [nGroup]=Convert2Ratios(Group)
sz=size(Group);
nGroup=zeros(sz);
for i=1:sz(1)
    base=sum(Group(i,:));
    nGroup(i,:)=Group(i,:)/base;
end
end
function [Map]=mnl_SelectRGBcolours(n)
if n==1
    Map=[1 0 0];
elseif n==2
    Map=[0 1 0;1 0 0];
elseif n==3
    Map=[0 0 1;0 1 0;1 0 0];
elseif n==4
    Map=[0 0 1;0 1 1;0 1 0;1 0 0];
elseif n==5
    Map=[0 0 1;0 1 1;0 1 0;1 1 0;1 0 0];
elseif n==6
    Map=[0.5 0 1;0 0 1;0 1 1;0 1 0;1 1 0;1 0 0];
elseif n==7
    Map=[0.5 0 1;0 0 1;0 1 1;0 1 0;1 1 0;1 0.5 0;1 0 0];
elseif n==8
    Map=[1 0 1;0.5 0 1;0 0 1;0 1 1;0 1 0;1 1 0;1 0.5 0;1 0 0];
end
end