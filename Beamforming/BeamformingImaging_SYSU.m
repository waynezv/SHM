%%% Wideband Beamforming for Imaging
%%% Author: Yuanwei Jin, University of Maryland Eastern Shore
%%% Created April 12, 2005 by Yuanwei Jin at Carnegie Mellon University
%%% Modified September 3, 2015 by Yuanwei Jin at UMES
%%% 
%%% Note: This code is to generate an image by wideband signaling
%%% Data is collected by experiments 
%%% Suggestions: This code was intended for radio frequency imaging in far field,
%%% however, the principle also applies to acoustic imaging, especially for
%%% Lamb waves on plates. The Lamb waves can be readily modelled by Hankel
%%% functions. 

clear all; 
close all
tic
%%
%% this line specifies the path, you can change 
% pp = 'C:\Users\yjin\Documents\Students\SYSU\';
addpath ./
%%
% load([pp 'frequency_range.mat']);
% load([pp 'Position_x.mat']);
% load([pp 'Position_y.mat']);
load frequency_range.mat
load Position_x.mat
load Position_y.mat
%% Running Case -8 
% load([pp 'Setup11.mat']);
% load([pp 'Setup13.mat']);
load Setup11.mat
load Setup13.mat
clutterID= [1 2 7  9 10 11 14 15 16 18 26 27 28 33 36 40 44];
targetID=[58];
Clutter=Setup11;
TargetClutter=Setup13;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define targets 
Numtarget=length(targetID);
Numscatter=length(clutterID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify the sensor arrays 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 10; %% N-antenna array 
%% frequency grid in GHz
frange=length(frequency_range);
fmin=min(frequency_range)/1e9;
fmax=max(frequency_range)/1e9;
fstep=(fmax-fmin)/(frange-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idA=1:N;
idB=1:N;
%% antenna locations in centi-meters --Y axis; X = 0;
TA=2.54*[48.6:-4:9]; % 2.54??
TB=-2.54*[10:4:48];
TA=TA(idA);
TB=TB(idB);
%%  the measurements
H3=reshape(Clutter, N,N, frange); 
H4=reshape(TargetClutter,N,N,frange);
K3=H3(idB,idA,:);
K4=H4(idB,idA,:);
Ksub=K4 - K3; %% this is the target response matrix, after clutter being subtracted
%% wavenumber computation 
Lambda=3e+8./frequency_range*100;  %!!! lambda in centimeters !!!
Wavenumber=2*pi./Lambda;
%% image size
Nx=200; 
Ny=160; 
delta=6/8; %% centi-meters
%%%%%%%%%%%%%%%%%%%%%%%%%%%
xindex_offset=200;
xarray=delta*[0:Nx-1]+xindex_offset;
yarray=delta*[-Ny/2:Ny/2-1];
%% compute Greens' functions (the Hankel functions)
%% k*|X-R| that computes the distance of array R and each pixel
%% of the grid of the image.
argA=zeros(Nx,Ny,N);
argB=zeros(Nx,Ny,N);
for m=1:N
    yA=yarray-TA(m);
    yB=yarray-TB(m);
    [YYA,XX]=meshgrid(yA,xarray);
    [YYB,XX]=meshgrid(yB,xarray);
    argA(:,:,m)=-sqrt(XX.^2+YYA.^2); % distance
    argB(:,:,m)=-sqrt(XX.^2+YYB.^2); 
end
%%%%% loop over frequencies %%%%%%%%%%%
Ind=[1:201];  % we use all the frequency samples under a wide band scenario 
              % Here a total of 201 frequency samples are used. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gA=zeros(Nx,Ny,N,length(Ind));
gB=gA;
%% compute the Hankel functions
for ee=1:length(Ind)
   gA(:,:,:,ee)=besselh(0,Wavenumber(Ind(ee))*argA);
   gB(:,:,:,ee)=besselh(0,Wavenumber(Ind(ee))*argB);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bkprop=zeros(Nx,Ny);
for ee=1:length(Ind)
   temp=zeros(Nx,Ny);
   for m=1:Nx
   for n=1:Ny
      vB=reshape(gB(m,n,:,ee),N,1);
      vA=reshape(gA(m,n,:,ee),N,1);
      vB=vB/sqrt(vB'*vB);  %% normalization
      vA=vA/sqrt(vA'*vA);  %% normalization 
      %% this is for beamforming 
      temp(m,n)=vB'*Ksub(:,:,Ind(ee))*vA;
   end;   
   end;
%% for conventional imaging
bkprop=bkprop+abs(temp).^2;
end;

%% Plotting 
bkpropdb=10*log10(bkprop);
MM=max(max(bkpropdb));
bkpropdb=bkpropdb-MM;
figure(10);
imagesc(xarray,yarray,bkpropdb.');
for ei=1:Numscatter
    H=text(Position_x(clutterID(ei)),Position_y(clutterID(ei)), num2str(clutterID(ei)));
   set(H,'FontSize',14);
end;
for ei=1:Numtarget
   H=text(Position_x(targetID(ei)), Position_y(targetID(ei)),'X');
   set(H,'FontSize',14);
   set(H,'color','white');
end;
xlabel('range [cm]');
ylabel('cross range [cm]');
colorbar;
title('Direct Subtraction Wideband Beamforming [dB]')

[dummy,ymax]=max(max((bkpropdb)));
[dummy,xmax]=max((bkpropdb(:,ymax)));
H=text(delta*(xmax-1)+min(xarray), delta*(ymax-1)+min(yarray),'o');
set(H,'FontSize',14);
set(H,'color','white');

processingtime=toc 



