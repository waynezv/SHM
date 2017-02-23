% This script is for acoustical beamforming imaging
clear all
% Add path 
addpath('../../ddmfp-tools')
addpath('../../dsp-tools')

%% SIMULATION CONFIGURATIONS        
% Define the location of sensors
% - PZT transducers on the left edge of a 1m*1m plate, 10 for sending and 10 for receiving
sendLoc = 2.25.*[50:-4:10]; % cm
sendLoc = [zeros(size(sendLoc,2),1) sendLoc.'];
receLoc = -2.25.*[10: 4:50];
receLoc = [zeros(size(receLoc,2),1) receLoc.'];
numSensor = size(sendLoc,1);
% Define the location of targets
targLoc = [50 0];
% - Put in a cell
location = {sendLoc, receLoc, targLoc};
% Define the signal parameters
% - Frequency of interest
freqRange = [10:1:1000]; % KHz
% - Sampled wavenumbers
k = [4:4:2000].'; % 1/m

% -Sample time
Fs = 10000000;
Q = 10000;
t = 1000*[1/Fs:1/Fs:Q/Fs]; 
% Define the wave propagation model
% - Model: X = D_rAD_k * V
% + X:The signal response of lamb waves, of (# of distances * # of
%     frequencies), lamb wave excited by 10 \mu s linear freqency modulated
%     chirp.
% ++ D_rAD_k is the sensing matrix, comprising of orthogonal measurement
%    matrices D_r and D_k, as well as basis matrix A
% + D_r:The distance matrix 
% + A:The basis of Fourier transform, e^j(kr), 
% + D_k:The wavenumber matrix
% + V: 

% - Define the distance matrix
% + distance between senders and receivers
d_1 = dist(location{1}, location{3}.'); % distances between senders and target
d_2 = dist(location{2}, location{3}.'); % distances between receivers and target
M = size(d_1, 1); % # of distances
mode = 'reflection'% select mode, 'reflection' for hitting the target and reflect back, distance computed as the addition of d_1 and d_2

switch mode
    case 'single' % 'single' for hitting or receiving only
        D_r = spdiags(1./sqrt(d_1), 0, M, M);
        % - Define the basis matrix
        A = exp(-1j*d_1*k.'); 
    case 'reflection'
        for path = 1:length(d_2)
            d(:,path) = d_1+repmat(d_2(path,:), [length(d_1), 1]);
        end
%         D_r = spdiags(1./sqrt(d), 0, M, M);
        % - Define the basis matrix
%         A = exp(-1j*d*k.'); 
end

% - Define the wavenumber matrix
N = size(k, 1); % # of wavenumbers
D_k = spdiags(1./sqrt(k), 0, N, N);
% - Introduce the dispersive matrix
load V_1M_s0.02_w1.5-0.5-0.25_t0.45.mat

% Compute signal response at the target from senders directly
switch mode
    case 'single'
        X =  D_r*A*D_k * V;
        x = ifft(X.');
    case 'reflection'
        for count = 1:size(d, 2)
            D_r = spdiags(1./sqrt(d(:,count)), 0, M, M);
            A = exp(-1j*d(:,count)*k.'); 
            X(count,:,:) = D_r*A*D_k * V;
        end
        for count = 1:size(d, 2)
            x(count,:,:) = ifft(X(count,:,:));
        end
end
% Save variables
s = struct('numSensor',{numSensor},'location',{location},'data',{x});
fprintf('saving variables...\n');
save('beamforming', '-struct', 's');
%% PLOT
% [X,Y] = meshgrid(t(1:size(x,3)), t(1:size(x,3)));
figure
for i = 1:numSensor
    subplot(6,2,i)
    plot(t(1:size(x,1)), x(:,i));
    title(['sensor channel' num2str(i)]);
end
suptitle('Signal responses')
% 
% figure,
% for i = 1:size(Rxm,2)
% plot(Rxm{i}, 'b+');
% hold on
% end
% for i = 1:size(Txm,2)
% plot(Txm{i}, 'ks');
% hold on
% end
% title('smdata2_line_ideal')