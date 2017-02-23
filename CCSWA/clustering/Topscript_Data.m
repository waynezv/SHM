clc;clear;
add_all_path();
% load('./Data/data_Aggregation.mat') %
% load('./Data/data_bridge.mat')
% load('./Data/data_Compound.mat')
% load('./Data/data_D31.mat')
% load('./Data/data_Flame.mat')
% load('./Data/data_Jain.mat')
% load('./Data/data_Pathbased.mat')
load('./Data/data_R15.mat')
% load('./Data/data_Spiral.mat')
% load('./Data/data_Target.mat')
% load('./Data/data_TwoDiamonds.mat')
point_matrix(:,1) = D(:,2);
point_matrix(:,2) = D(:,1);


cluster_num = 15;


%% Data Generation
% img = double(imread('Gaussian', 'bmp'));
% [nrow ncol channel] = size(img);
% total_point_num = sum(sum(img(:,:,1)<128));
% snr = -26;
% ratio = 400;
% point_matrix = zeros(total_point_num*2, 2);
% count_point = 1;
% generator_num = 0;
% for rr = 1 : nrow
%     for cc = 1 : ncol
%         if img(rr, cc, 1) < 128
%             generator_num = generator_num + 1;
%             if generator_num == 2 || generator_num == 4
%                 ratio = 500;
%                 snr = - 28;
%             else
%                 ratio = 250;
%                 snr = - 23;
%             end
%             for i = 1 : ratio
%                 point_matrix(count_point, 1:2) = [rr+awgn(0, snr) cc+awgn(0, snr)];
%                 count_point = count_point + 1;
%             end
%         end
%     end
% end

point_num = size(point_matrix, 1);
Dist_Mat= pdist2(point_matrix, point_matrix, 'euclidean');
% Graph_Edge = Build_RAG_Data(point_matrix);
% [~, Trans_Dist_Mat, ~] = MST(Graph_Edge);

Trans_Dist_Mat1 = MST_Data(Dist_Mat);
Trans_Dist_Mat1 = max_var_anal2(Trans_Dist_Mat1, cluster_num);

tree_num = 20;
% Trans_Dist_Mat2 = MSF_Data_Pert(Dist_Mat, tree_num, 2);
% Trans_Dist_Mat2 = mean(Trans_Dist_Mat_Set,3);
% for i = 1:point_num
%     Trans_Dist_Mat2(i,:) = Trans_Dist_Mat2(i,:)./sum(Trans_Dist_Mat2(i,:).^2);
% end

Trans_Dist_Mat2 = MSF_Data2(Dist_Mat, 5, 'max');
% Trans_Dist_Mat2 = Trans_Dist_Mat2*normrnd(0,1,[size(Trans_Dist_Mat2,1),15]);
Trans_Dist_Mat2 = max_var_anal2(Trans_Dist_Mat2, cluster_num);


% Trans_Dist_Mat = zeros(point_num, point_num);
% for i = 1:point_num
%     for j = 1:point_num
%         Trans_Dist_Mat(i, j) = norm(point_matrix(i, :)-point_matrix(j, :));
%         Trans_Dist_Mat(j, i) = norm(point_matrix(i, :)-point_matrix(j, :));
%     end
% end

% IDX1 = kmeans(Trans_Dist_Mat1,cluster_num, 'start', 'cluster', 'emptyaction', 'singleton');
% IDX2 = kmeans(Trans_Dist_Mat2,cluster_num, 'start', 'cluster', 'emptyaction', 'singleton');
% [IDX2, ~] = kmedoids(Trans_Dist_Mat2, cluster_num);
% IDX2 = IDX2';
IDX1 = kmeans3(Trans_Dist_Mat1,cluster_num,20);
IDX2 = kmeans3(Trans_Dist_Mat2,cluster_num,20);

figure(1)
for i = 1:cluster_num
    idx = IDX1==i;
    plot(point_matrix(idx, 1), point_matrix(idx, 2), '.', 'color', rand(1,3));hold on
end

figure(2)
for i = 1:cluster_num
    idx = IDX2==i;
    plot(point_matrix(idx, 1), point_matrix(idx, 2), '.', 'color', rand(1,3));hold on
end

% figure(1)
% for i = 1:size(point_matrix, 1)
%     if IDX1(i, 1) == 1
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'red');hold on
%     elseif IDX1(i, 1) == 2
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'green');hold on
%     elseif IDX1(i, 1) == 3
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'blue');hold on
%     elseif IDX1(i, 1) == 4
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'yellow');hold on
%     elseif IDX1(i, 1) == 5
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'cyan');hold on
%     elseif IDX1(i, 1) == 6
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'magenta');hold on
%     else
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'black');hold on
%     end
% end
% hold off
% 
% figure(2)
% for i = 1:size(point_matrix, 1)
%     if IDX2(i, 1) == 1
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'red');hold on
%     elseif IDX2(i, 1) == 2
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'green');hold on
%     elseif IDX2(i, 1) == 3
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'blue');hold on
%     elseif IDX2(i, 1) == 4
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'yellow');hold on
%     elseif IDX2(i, 1) == 5
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'cyan');hold on
%     elseif IDX2(i, 1) == 6
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'magenta');hold on
%     else
%         plot(point_matrix(i, 1), point_matrix(i, 2), '.', 'color', 'black');hold on
%     end
% end
% hold off
