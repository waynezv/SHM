clear; clc;
img = double(imread('Gaussian', 'bmp'));
[nrow ncol channel] = size(img);
% total_point_num = sum(sum(img(:,:,1)<128));
% snr = -26;
% ratio = 400;
% point_matrix = zeros(total_point_num*2, 2);
count_point = 1;
generator_num = 0;
for rr = 1 : nrow
    rr
    for cc = 1 : ncol
        if img(rr, cc, 1) < 128
            generator_num = generator_num + 1;
            if generator_num == 2 || generator_num == 4
                ratio = 500;
                snr = - 28;
            else
                ratio = 250;
                snr = - 23;
            end
            for i = 1 : ratio
                point_matrix(count_point, 1:2) = [rr+awgn(0, snr) cc+awgn(0, snr)];
                count_point = count_point + 1;
            end
        end
    end
end
total_point_num = count_point-1;
for i = 1 : total_point_num
    plot(point_matrix(i, 2), nrow - point_matrix(i, 1), '.'); hold on;
end
hold off;