clc; clear;
add_all_path();

%% Loading Data
% target_data_matfile='model.data.mat';
% model_target_name_matfile='target_speaker_models.csv';

load('./Data/data_dev_ivector.mat');
% load('./Data/data_dev_nist.mat')
% [data_dev,mu_dev,std_dev] = zscore_normalize(data_dev);% Z-score Normalization
load('./Data/label_dev_ivector.mat');
% load('./Data/label_dev_nist.mat');

% data_dev = data_dev(:,1:6000);
% label_dev_gt = label_dev_gt(1:6000);

%% Iterative Unsupervised Training
data_dev2 = data_dev;
for i=1:size(data_dev2,2)
    data_dev2(:,i)=data_dev2(:,i)./norm(data_dev2(:,i),2);
end
N = size(data_dev2,2);

Total_Clust_Num = 4958;
label_dev = kmeans(data_dev2',Total_Clust_Num,'distance','cosine');
score_purity_kms = cluster_purity(label_dev_gt,label_dev)
save('./intermediate_result/label_ivector_kms.mat', 'label_dev');