clc; clear;
add_all_path();

%% Loading Data
% target_data_matfile='model.data.mat';
% model_target_name_matfile='target_speaker_models.csv';

% load('./Data/data_dev_ivector.mat');
load('./Data/data_dev_nist.mat')
% [data_dev,mu_dev,std_dev] = zscore_normalize(data_dev);% Z-score Normalization
% load('./Data/label_dev_ivector.mat');
load('./Data/label_dev_nist.mat');

% data_dev = data_dev(:,1:6000);
% label_dev_gt = label_dev_gt(1:6000);

%% Iterative Unsupervised Training
data_dev2 = data_dev;
for i=1:size(data_dev2,2)
    data_dev2(:,i)=data_dev2(:,i)./norm(data_dev2(:,i),2);
end
N = size(data_dev2,2);

Total_Clust_Num = 1738;
K = 1;
idx_kms = kmeans(data_dev2',K,'distance','cosine','start','cluster');
idx_sub = cell(K,1);
data_dev_sub = cell(K,1);
cluster_num = cell(K,1);
count_total = 0;
for k = 1:K
    idx_sub{k,1} = idx_kms==k;
    data_dev_sub{k,1} = data_dev(:,idx_sub{k,1});
    cluster_num{k,1} = round(size(data_dev_sub{k,1},2)/(N/Total_Clust_Num));
    count_total = count_total+cluster_num{k,1};
end
if count_total ~= Total_Clust_Num %Make sure the total number of clusters stays the same with the original one. This is needed for NMI.
    num_diff = Total_Clust_Num-count_total;
    cluster_num{K,1} = cluster_num{K,1}+num_diff;
end

display('Computing pairwise similarity matrix...')
tic
dist_mat = cell(K,1);
similarity = cell(K,1);
for k = 1:K
    % dist_mat{k,1} = 0.5.*(1-data_dev2(:,idx_sub{k,1})'*data_dev2(:,idx_sub{k,1}));
    dist_mat{k,1} = pdist2(data_dev2(:,idx_sub{k,1})', data_dev2(:,idx_sub{k,1})');
    similarity{k,1} = 1-dist_mat{k,1};
end
toc
display(['Computing pairwise dist_mat matrix takes ' num2str(toc) ' seconds']);

label_dev_trans1 = zeros(N,1);
label_dev_trans2 = zeros(N,1);
label_dev_RF_MultiPath1 = zeros(N,1);
label_dev_RF_MultiPath2 = zeros(N,1);
label_dev_RF_Pert1 = zeros(N,1);
label_dev_RF_Pert2 = zeros(N,1);
label_dev_Spectral = zeros(N,1);
label_dev_Ncut = zeros(N,1);
label_dev_AP = zeros(N,1);
label_dev_AHC = zeros(N,1);
label_dev_Kmeans = zeros(N,1);
pointer = zeros(11, 1);

for k = 1:K
    k
    %% MST Transitive Distance (No Dim Reduction)
    trans_dist_mat = MST_Data(dist_mat{k,1});
    trans_dist_mat = trans_dist_mat.*dist_mat{k,1};
    %     tic;
    %     label_dev_tmp = kmeans(trans_dist_mat, cluster_num{k,1}, 'emptyaction', 'singleton', 'distance', 'cosine');
    %     display(['MST Trans Kmeans (No Dim Reduction) takes ' num2str(toc) ' seconds']);
    %     label_dev_tmp = unify_label(label_dev_tmp) + pointer(1,1);
    %     label_dev_trans1(idx_sub{k,1},1) = label_dev_tmp;
    %     pointer(1,1) = pointer(1,1) + length(unique(label_dev_tmp));
    
    %% MST Transitive Distance (Dim Reduction)
    %[trans_dist_mat] = max_var_anal(trans_dist_mat, min(1500,size(trans_dist_mat,1)));
    tic
    trans_dist_mat = max_var_anal2(trans_dist_mat, 1400);
    display(['MST Trans SVD takes (Dim Reduction) ' num2str(toc) ' seconds']);
    tic;
    purity_set = zeros(11,1);
    counter = 0;
    for column_num = [100 200 300 400 500 1400]
        counter = counter+1;
        column_num
        label_dev_tmp = kmeans(trans_dist_mat(:,1:column_num), cluster_num{k,1}, 'emptyaction', 'singleton', 'distance', 'cosine');
        score_purity_trans2 = cluster_purity(label_dev_gt,label_dev_tmp)
        purity_set(counter,1) = score_purity_trans2;
        label_dev = label_dev_tmp;
        save(['./intermediate_result/label_dev_nist_lctd2_column_num_' num2str(column_num) '.mat'], 'label_dev');
    end
    % save('./purity_dev_nist_lctd2.mat', 'purity_set');
    % label_dev_tmp = kmeans4(trans_dist_mat, cluster_num{k,1}, 10);
    display(['MST Trans Kmeans takes (Dim Reduction) ' num2str(toc) ' seconds']);
    label_dev_tmp = unify_label(label_dev_tmp) + pointer(2,1);
    label_dev_trans2(idx_sub{k,1},1) = label_dev_tmp;
    pointer(2,1) = pointer(2,1) + length(unique(label_dev_tmp));
    
    %     %% Random Forest Transitive Distance (Multipath, No Dim Reduction)
    %     trans_dist_mat = MSF_Data2(dist_mat{k,1}, 5, 'max');
    %     tic;
    %     label_dev_tmp = kmeans(trans_dist_mat, cluster_num{k,1}, 'emptyaction', 'singleton', 'distance', 'cosine');
    %     display(['MSF Multipath Kmeans (No Dim Reduction) takes ' num2str(toc) ' seconds']);
    %     label_dev_tmp = unify_label(label_dev_tmp) + pointer(3,1);
    %     label_dev_RF_MultiPath1(idx_sub{k,1},1) = label_dev_tmp;
    %     pointer(3,1) = pointer(3,1) + length(unique(label_dev_tmp));
    %
    %     %% Random Forest Transitive Distance (Multipath, Dim Reduction)
    %     [trans_dist_mat] = max_var_anal(trans_dist_mat, min(1500,size(trans_dist_mat,1)));
    %     tic
    %     trans_dist_mat = max_var_anal2(trans_dist_mat, cluster_num{k,1});
    %     display(['MSF Multipath SVD takes (Dim Reduction) ' num2str(toc) ' seconds']);
    %     tic;
    %     label_dev_tmp = kmeans(trans_dist_mat, cluster_num{k,1}, 'emptyaction', 'singleton', 'distance', 'cosine');
    %     display(['MSF Multipath Kmeans (Dim Reduction) takes ' num2str(toc) ' seconds']);
    %     label_dev_tmp = unify_label(label_dev_tmp) + pointer(4,1);
    %     label_dev_RF_MultiPath2(idx_sub{k,1},1) = label_dev_tmp;
    %     pointer(4,1) = pointer(4,1) + length(unique(label_dev_tmp));
    %
    %     %% Random Forest Transitive Distance (Perturbation, No Dim Reduction)
    %     trans_dist_mat = MSF_Data_Pert(dist_mat{k,1}, 20, 0.015);
    %     tic;
    %     label_dev_tmp = kmeans(trans_dist_mat, cluster_num{k,1}, 'emptyaction', 'singleton', 'distance', 'cosine');
    %     display(['MSF Pert Kmeans (No Dim Reduction) takes ' num2str(toc) ' seconds']);
    %     label_dev_tmp = unify_label(label_dev_tmp) + pointer(5,1);
    %     label_dev_RF_Pert1(idx_sub{k,1},1) = label_dev_tmp;
    %     pointer(5,1) = pointer(5,1) + length(unique(label_dev_tmp));
    %
    %     %% Random Forest Transitive Distance (Perturbation, Dim Reduction)
    % %     [trans_dist_mat] = max_var_anal(trans_dist_mat, min(1500,size(trans_dist_mat,1)));
    %     tic
    %     trans_dist_mat = max_var_anal2(trans_dist_mat, cluster_num{k,1});
    %     display(['MSF Pert SVD (Dim Reduction) takes ' num2str(toc) ' seconds']);
    %     tic;
    % %     label_dev_tmp = kmeans(trans_dist_mat, cluster_num{k,1}, 'emptyaction', 'singleton', 'distance', 'cosine');
    %     label_dev_tmp = kmeans4(trans_dist_mat, cluster_num{k,1}, 10);
    %     display(['MSF Pert Kmeans (Dim Reduction) takes ' num2str(toc) ' seconds']);
    %     label_dev_tmp = unify_label(label_dev_tmp) + pointer(6,1);
    %     label_dev_RF_Pert2(idx_sub{k,1},1) = label_dev_tmp;
    %     pointer(6,1) = pointer(6,1) + length(unique(label_dev_tmp));
    %
    %     %% Spectral Clustering
    %     tic;
    %     label_dev_tmp = spectral_clust2(similarity{k,1}, cluster_num{k,1}, cluster_num{k,1});
    %     label_dev_tmp = unify_label(label_dev_tmp) + pointer(7,1);
    %     label_dev_Spectral(idx_sub{k,1},1) = label_dev_tmp;
    %     pointer(7,1) = pointer(7,1) + length(unique(label_dev_tmp));
    %     display(['Spectral clustering takes ' num2str(toc) ' seconds']);
    %
    %     %% Normalized Cuts
    %     tic;
    %     label_dev_tmp = ncut_clust(similarity{k,1}, cluster_num{k,1});
    %     label_dev_tmp = unify_label(label_dev_tmp) + pointer(8,1);
    %     label_dev_Ncut(idx_sub{k,1},1) = label_dev_tmp;
    %     pointer(8,1) = pointer(8,1) + length(unique(label_dev_tmp));
    %     display(['Ncut takes ' num2str(toc) ' seconds']);
    
    %     %% Affinity Propagation
    %     tic;
    %     [label_dev_tmp,~,~,~] = apcluster(similarity{k,1},median(similarity{k,1}),'maxits',200);
    %     label_dev_tmp = unify_label(label_dev_tmp) + pointer(9,1);
    %     label_dev_AP(idx_sub{k,1},1) = label_dev_tmp;
    %     pointer(9,1) = pointer(9,1) + length(unique(label_dev_tmp));
    %     display(['AP takes ' num2str(toc) ' seconds']);
    
    %     %% Single Linkage (Agglomerative Hierarchical Clustering)
    %     tic;
    %     label_dev_tmp = AHC(dist_mat{k,1}, cluster_num{k,1});
    %     label_dev_tmp = unify_label(label_dev_tmp)+pointer(10,1);
    %     label_dev_AHC(idx_sub{k,1},1) = label_dev_tmp;
    %     pointer(10,1) = pointer(10,1) + length(unique(label_dev_tmp));
    %     display(['AHC takes ' num2str(toc) ' seconds']);
    %
    %     %% Kmeans
    %     tic;
    %     label_dev_tmp = kmeans(data_dev2(:,idx_sub{k,1})', cluster_num{k,1}, 'emptyaction', 'singleton', 'distance', 'cosine');
    %     label_dev_tmp = unify_label(label_dev_tmp) + pointer(11,1);
    %     label_dev_Kmeans(idx_sub{k,1},1) = label_dev_tmp;
    %     pointer(11,1) = pointer(11,1) + length(unique(label_dev_tmp));
    %     display(['Kmeans takes ' num2str(toc) ' seconds']);
end

% %% MST Transitive Distance (No Dim Reduction)
% score_purity_trans1 = cluster_purity(label_dev_gt,label_dev_trans1)
% 
% %% MST Transitive Distance (Dim Reduction)
% score_purity_trans2 = cluster_purity(label_dev_gt,label_dev_trans2)

% %% Random Forest Transitive Distance (Multipath, No Dim Reduction)
% score_purity_multipath1 = cluster_purity(label_dev_gt,label_dev_RF_MultiPath1)
%
% %% Random Forest Transitive Distance (Multipath, Dim Reduction)
% score_purity_multipath2 = cluster_purity(label_dev_gt,label_dev_RF_MultiPath2)
%
% %% Random Forest Transitive Distance (Perturbation, No Dim Reduction)
% % score_purity_pert1 = cluster_purity(label_dev_gt,label_dev_RF_Pert1)
%
% %% Random Forest Transitive Distance (Perturbation, Dim Reduction)
% score_purity_pert2 = cluster_purity(label_dev_gt,label_dev_RF_Pert2)
%
% %% Spectral Clustering
% score_purity_spectral = cluster_purity(label_dev_gt,label_dev_Spectral)
%
% %% Normalized Cuts
% score_purity_ncut = cluster_purity(label_dev_gt,label_dev_Ncut)
%
% %% Affinity Propagation
% % score_purity_ap = cluster_purity(label_dev_gt,label_dev_AP)
%
% %% Single Linkage (Agglomerative Hierarchical Clustering)
% score_purity_ahc = cluster_purity(label_dev_gt,label_dev_AHC)
%
% %% Kmeans
% score_purity_kms = cluster_purity(label_dev_gt,label_dev_Kmeans)



% max_iter = 1;
% beta = 0.2;
%
% EER_Set = zeros(max_iter, 2);
% DCF_Set = zeros(max_iter, 2);
% for j=1:max_iter
%     %% Construct Training Data
%     spklabel_gplda=[];
%     data_dev_train = [];
%     for k=1:K
%         spklabel_gplda_idx=[];
%         spklabel=unique(label_dev{k,1});
%         spklabelidx=length(spklabel);
%         for i=1:spklabelidx
%             idx=find(label_dev{k,1}==spklabel(i));
%             spklabel_gplda_idx=[spklabel_gplda_idx;idx];
%             spklabel_gplda=[spklabel_gplda;length(idx)];
%         end
%         data_dev_train=[data_dev_train data_dev_sub{k,1}(:,spklabel_gplda_idx)];
%     end
%
%     %% Train PLDA Model
%     [Scores,m,W,Phi,Sigma] = Gaussian_PLDA(data_dev_train,spklabel_gplda,data_dev,data_dev,600,1);
%
%     %% Construct Target Set
%     data_target2 = bsxfun(@minus,data_target,m);
%     data_test2 = bsxfun(@minus,data_test,m);
%     data_target2 = W*data_target2;
%     data_test2 = W*data_test2;
%     m_nrm = sqrt(sum(data_target2 .* data_target2));
%     t_nrm = sqrt(sum(data_test2 .* data_test2));
%     data_target2 = bsxfun(@rdivide,data_target2,m_nrm);
%     data_test2 = bsxfun(@rdivide,data_test2,t_nrm);
%
%     data_model=[];
%     for i=1:idx_model
%         idx=find(label_target==i);
%         the_target_data=data_target2(:,idx);
%         assert(length(idx)==5);
%         the_target_data=reshape(the_target_data,size(the_target_data,2)*size(the_target_data,1),1);
%         data_model=[data_model the_target_data];
%         clear idx;
%         clear the_target_data;
%     end
%
%     %% Compute Final Score
%     thescores=PLDA_5target1test_score(Phi,Sigma,data_model,data_test2);
%
%     t1=find(ivec14_sre_trial_key(:,2)==1);%prog
%     t2=find(ivec14_sre_trial_key(t1,1)==0);%nontarget
%     t3=find(ivec14_sre_trial_key(t1,1)==1);%target
%
%     prog_true=t1(t3);
%     prog_false=t1(t2);
%     clear t1;
%     clear t2;
%     clear t3;
%
%     t1=find(ivec14_sre_trial_key(:,2)==2);%eval
%     t2=find(ivec14_sre_trial_key(t1,1)==0);%nontarget
%     t3=find(ivec14_sre_trial_key(t1,1)==1);%target
%
%     eval_true=t1(t3);
%     eval_false=t1(t2);
%
%     clear t1;
%     clear t2;
%     clear t3;
%
%     finalscores=reshape(thescores',size(thescores,2),size(thescores,1));
%     [EER_prog, DCF_opt_prog]=getresult2(finalscores(prog_true),finalscores(prog_false),1);
%     [EER_eval, DCF_opt_eval]=getresult2(finalscores(eval_true),finalscores(eval_false),2);
%     EER_prog
%     EER_eval
%     DCF_opt_prog
%     DCF_opt_eval
%
%     EER_Set(j,:) = [EER_prog EER_eval];
%     DCF_Set(j,:) = [DCF_opt_prog DCF_opt_eval];
%     save('./Result_transdist_sub_znorm.mat', 'EER_Set', 'DCF_Set');
%
%     %% Incremental Update with Iteration
%     if j==max_iter
%         break
%     else
%         t1=min(min(Scores));
%         t2=max(max(Scores));
%         Scores=(Scores-t1)./(t2-t1);% Normalize Score
%
%         for k = 1:K
%             dist_mat{k,1}=(1-Scores(idx_sub{k,1},idx_sub{k,1})).*beta+dist_mat{k,1}.*(1-beta);
%         end
%
%         display('Using Subclass Spectral clustering to perform clustering...');
%         for k = 1:K
%             k
%             label_dev{k,1} = spectral_clust2(dist_mat{k,1}, 600, cluster_num{k,1});
%         end
%         display(['Subclass Spectral clustering takes ' num2str(ttt) ' seconds']);
%     end
% end