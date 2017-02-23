% ------ Wideband Lamb Wave Beamforming for Imaging ------ %
% % 
% %
% -
% -

clear all, close all, clc
tic

%% Data Settings
addpath ./;
load beamforming.mat
% Pre-configuration
% - wavenumber
k = [4:4:2000].'; % 1/m
% - frequency
freqRange = [10:1:1000]; % KHz
% - number of sensors
N = numSensor;
% -Sample time
Fs = 10000000;
Q = 10000;
t = 1000*[1/Fs:1/Fs:Q/Fs]; 
% - distance
d_s2t = dist(location{1}, location{3}.'); % distances between senders and target
d_t2r = dist(location{2}, location{3}.'); % distances between receivers and target
d = d_s2t + d_t2r;
%%  The measurementss
% - signal response at the receivers
x = data;
for i = 1:size(x,1)
    X(i,:,:) = fft(x(i,:,:)); % 11*991
end
%% Learn the dispersive properties of plate
% + options for SWA
 method = 'omp';          % Optimization algorithm used 
 tau    = 3;              % BP regularization parameter or OMP sparsity
% UNCOMMENT THE NEXT TWO LINES TO USE BASIS PURSUIT DENOISING
% method = 'bp';           % Optimization algorithm used
% tau    = 0.45;         % BP regularization parameter or OMP sparsity

% method='mli';
% tau    = 0.45;

optMatrix    = false;    % Solve the basis pursuit for all frequencies at once
optNomalized = false;   % Use distortionless constraint
iterPlot     = false;    % Plot result at each frequency
% + SWA
% [V_s2t, Phi] = swa_mli_v2_ori( k, d_s2t, X(, tau,'optMatrix', optMatrix, 'optNomalized', optNomalized, 'plot', iterPlot, 'method', method );
% [V_t2r, ~] = swa_mli_v2_ori( k, d_t2r, X, tau,'optMatrix', optMatrix, 'optNomalized', optNomalized, 'plot', iterPlot, 'method', method );
for path = 1:length(d_t2r)
    d(:,path) = d_s2t+repmat(d_t2r(path,:), [length(d_s2t), 1]);
end
d = d(:);

for i = 1:size(X,1)
    for j = 1:size(X,2)
        X2(i*j,:) = X(i,j,:);
    end
end

[V, ~] = swa_mli_v2_ori( k, d, X2.', tau,'optMatrix', optMatrix, 'optNomalized', optNomalized, 'plot', iterPlot, 'method', method );
% + SWS
x_synthesis = 2*real(sws_mli_v1( k, d, V ));
%% Image Settings
Nx=200; 
Ny=160; 
delta=6/8; %% centi-meters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
xindex_offset=200;
xarray=delta*[0:Nx-1]+xindex_offset;
yarray=delta*[-Ny/2:Ny/2-1];
% - grid
% k*|X-R| that computes the distance of array R and each pixel
% of the grid of the image.
argA=zeros(Nx,Ny,N);
argB=zeros(Nx,Ny,N);
TA = cell2mat(location(1));
TA = TA(:,2).';
TB = cell2mat(location(2));
TB = TB(:,2).';
for m=1:N
    yA=yarray-TA(m);
    yB=yarray-TB(m);
    [YYA,XX]=meshgrid(yA,xarray);
    [YYB,XX]=meshgrid(yB,xarray);
    argA(:,:,m)=-sqrt(XX.^2+YYA.^2);
    argB(:,:,m)=-sqrt(XX.^2+YYB.^2); 
end

%% Compute Greens' functions (the Hankel functions)
%%%%% loop over frequencies %%%%%%%%%%%
Ind=[1:size(X,3)];  % we use all the frequency samples under a wide band scenario 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gA=zeros(Nx,Ny,N,length(Ind));
gB=gA;
% Wavenumber = zeros(size(V));
for j = 1:length(Ind)
    Wavenumber(:,j) = k(find(abs(V(:,j))~=0));
end

for ee=1:length(Ind)
   gA(:,:,:,ee)=besselh(0,Wavenumber(Ind(ee))*argA);
   gB(:,:,:,ee)=besselh(0,Wavenumber(Ind(ee))*argB);
end
%% Formulate the beamforming output
bkprop=zeros(Nx,Ny);
for ee=1:length(Ind)
   temp=zeros(Nx,Ny);
   for m=1:Nx
   for n=1:Ny
      vB=reshape(gB(m,n,:,ee),N,1);
      vA=reshape(gA(m,n,:,ee),N,1);
      vB=vB/sqrt(vB'*vB);  %% normalization
      vA=vA/sqrt(vA'*vA);  %% normalization 
      % this is for beamforming 
      temp(m,n)=vB'*X(:,:,Ind(ee))*vA;
   end   
   end
% for conventional imaging
bkprop=bkprop+abs(temp).^2;
end

%% Plot
bkpropdb=10*log10(bkprop);
MM=max(max(bkpropdb));
bkpropdb=bkpropdb-MM;
figure(10);
imagesc(xarray,yarray,bkpropdb.');

figure
imagesc(abs(V), [0 0.5]);

figure
plot(x_synthesis(:,2));

processingtime=toc 