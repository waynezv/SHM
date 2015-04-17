[BW,thresh,gv,gh] = edge(abs(Vn),'sobel',0.02);

X_axis = [1:1:500];
Y_axis = [1:1:991];

Y_axis2 = [1:1:1000];

% C(:,:,1)=R; C(:,:,2)=G;C(:,:,3)=B; 
% C = zeros(size(500, 1000));
h = mesh(Y_axis, X_axis, abs(Vn));
set(h,'EdgeColor','k','FaceColor','w','MarkerEdgecolor','b','MarkerFacecolor','c');

subplot(131), mesh(Y_axis, X_axis, abs(Vn));
subplot(132), mesh(Y_axis, X_axis, gv);
subplot(133), mesh(Y_axis, X_axis, gh);

Vnsub = abs(Vn)-gv./abs(gv).*mean(abs(Vn(:)));

figure(2), mesh(Y_axis, X_axis, abs(Vnsub));

% plot clusters
orimap=abs(Vn)>0.02;
feat=[];
for i=1:size(orimap,1)
    for j=1:size(orimap,2)
        if(orimap(i,j)>0)
            feat=[feat;[i,j]];
        end
    end
end

feat2=[feat(:,1) feat(:,2)-feat(:,1)];
x=feat2';
% [a, model, L] = vbgm_wz_1(x, 10);
% label_name=unique(a);
label=vbgm_wz_1(x,10); 
spread(flipud(x),label);

surf(Y_axis, X_axis, abs(Vn));

mesh(Y_axis, X_axis, abs(Vn));

contour(Y_axis, X_axis, abs(Vn));

% Gaussian blur and unsharpen mask
[c, h] = contour(abs(Vn > mean(Vn(:))));
saveas(h, '');

% if in a loop
saveas(gcf,['PATH','FILENAME',num2str(i),'.jpg']); % i indicates the number of ith loop

% or 
test2 = getframe(gcf);
imwrite(test2.cdata,'test2.jpg');
imwrite(A, FILENAME, FMT);

% Blur
A = imread('test1.jpg');
gFil = fspecial('gaussian',[5 5],2);
gBlur = imfilter(A, gFil); 
imshow(gBlur);
saveas(gcf, 'test1Blr.jpg');
% sharpen
I = imread('test1Blr.jpg');
Ishrp = imsharpen(I, 'Radius', 2, 'Amount', 1);
imshow(Ishrp);
% substitution of image unsharpen
I = imread('test1Blr.jpg');
gUnshrp = padarray(2,[2 2]) - fspecial('gaussian' ,[5 5],2); % create unsharp mask
Ishrp2 = imfilter(I, gUnshrp);
imshow(Ishrp2);
% substitution of implementing gaussian blur
    I1 = fft2 (I);
    I2 = fftshift (I1);
    [w, h] = size (I);
    GF = zeros ([w, h]);
    % This value determines the 'blur' amount:
    % The higher it goes, the less the blur is;
    % The lower it goes, the more blur there is.
    Sig = 5;
    for x = 1 : w
        for y = 1 : h
            GF (x, y) = exp (-((x-w/2)^2 + (y-h/2)^2) / (2 * Sig^2));
        end
    end
    I3 = I2 .* GF;
    I4 = ifft2 (fftshift (I3));

% functions for plot
scatter(X, Y, S, C); % 2D or 3D

imshow();

imagesc(X, clims); % 2D or 3D

contour();

mesh(); % with or without meshgrid

surf();