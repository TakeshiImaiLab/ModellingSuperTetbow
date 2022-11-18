function [EuD_all]=mnl_GroupColourEuclidean_Simplified(data)
%% Help Section
%
% Inputs
% data - data points as vector (cells*ndim)
%
% Outputs
% EuD_all - Euclidean Distance between all points
%
% Written by Marcus Leiwe, Kyushu University, 2018. If used please cite
% Sakaguchi et al 2018

%% Basic Info
sz=size(data);
EuD_all=nan((sz(1)-1)*(sz(1)/2),1);
n=1;
%% Calculate Euclidean Distances
for i=1:(sz(1)-1) %For Each Row
    for j=i+1:sz(1)
        val=MNL_EuclideanDistance(data(i,:),data(j,:));
        EuD_all(n,1)=val;
        n=n+1;
    end
end
end

function [EuD]=MNL_EuclideanDistance(X,Y)
EuD=sqrt(sum((X-Y).^2));
end