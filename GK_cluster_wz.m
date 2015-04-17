% GK clustering feat3

clear all
close all
addpath('.\FUZZCLUST');

load Vn_polyfit.mat;
colors={'r.' 'gx' 'b+' 'ys' 'md' 'cv' 'k.' 'r*' 'g*' 'b*' 'y*' 'm*' 'c*' 'k*' };

%the data
% data.X=nDexample(5,250,2,1);
data.X = Vn;
%normalization
data=clust_normalize(data,'range');

%parameters
param.c=6;
param.m=2;
param.e=1e-3;
param.ro=ones(1,param.c);
param.beta = 1e15;
param.gamma = 0;
param.val=1;

%Gustafson Kessel-clustering
result = GKclust(data,param);

%validation
result = validity(result,data,param);
plot(data.X(:,1),data.X(:,2),'b.',result.cluster.v(:,1),result.cluster.v(:,2),'ro');
hold on
plot(result.cluster.v(:,1),result.cluster.v(:,2),'ro');

%evaluation
new.X=data.X;
eval = clusteval(new,result,param);
result.validity