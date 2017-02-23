% -------------------------------------------------------------------------
% EXAMPLE_SWA_SWS
% Written by: Joel Harley
% Last updated: May 5, 2014
% -------------------------------------------------------------------------
%
% In this script, we show an example of applying sparse wavenumber
% analysis (SWA) to simulation data corrupted by multipath interference.  
%
% The file "simdata_ideal" and "simdata_multipath" contains the ideal 
% (no multipath) or multipath "measured" simulation data. This data is used
% as the input to sparse wavenumber analysis (SWA). 
%
% The file "simdata" contains the ideal (no multipath) simulation data. 
% This data is used to evaluate our sparse wavenumber denoising results -- 
% i.e., test if the denoised data matches the true ideal deta. 
%
% The file "simdata_predict" contains a set of ideal responses from 
% random locations in the medium. This data is used to evaluate our sparse 
% wavenumber prediction results -- i.e., can we accurately predict other
% "measured" signals from the medium. 
%
% NOTE: This script is signifigantly slower than "example_simdata" since
% more data is necessary to overcome the multipath corruption. 
%
% 
% The script uses two particular functions: SWA and SWS.
%
% SWA performs sparse wavenumber analysis on the data 'xm' at frequency
% samples 'fn' and over a wavenumber domain 'k'. Sparse wavenumber analysis
% can be performed using one of two methods currently: basis pursuit 
% denoising ('bp') or orthogonal matching pursuit ('omp'). 
% 
% The variables 'Rxm' and 'Txm' represent the looks of the recievers and
% transmitters, respectively. Each cell in 'Rxm' and 'Txm' represent a 
% different experimental trial. The variable 'tau' represents the 
% regulerization parameter used by basis pursuit denoising OR the number 
% of sparse components for orthogonal matching pursuit. The output 'V' 
% is a matrix of the data's frequency-wavenumber representation. 
%
% SWS performs sparse wavenumber synthesis given a frequency-wavenumber
% representation 'V' over a wavenumber domain 'k'. The variables 'Rxm' and 
% 'Txm' or 'Rxp' and 'Txp' represent the virtual locations sensors or
% predicting each response. 
%
% For more information on the concepts behing this code, see: 
%   Sparse Recovery of the Multimodal and Dispersive  
%   Characteristics of Lamb Waves
%   J.B. Harley, J.M.F. Moura
%   Journal of the Acoustical Society of America
%   vol. 133, no. 5, pp. 2732-2745, May 2013
%
% -------------------------------------------------------------------------

close all;
% LOAD DATA (UNCOMMENT ONE OF THE SETUPS BELOW)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% DATA: RANDOMLY SEPARATED SENSORS WITH MULTIPATH (max M = 240)
addpath ../swa_data
addpath ../data_mat_t/
addpath('../cvx-w64');
addpath ./TVnorm
%   Note: May want to use around M = 200 measurements for this
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load simdata1_multipath; % Multipath data 
load simdata1;           % Ideal data 
load simdata2_predict;   % Data to predict (same as ideal data right now)
%  method = 'omp';          % Optimization algorithm used 
%  tau    = 3;              % BP regularization parameter or OMP sparsity
% UNCOMMENT THE NEXT TWO LINES TO USE BASIS PURSUIT DENOISING
% method = 'bp';           % Optimization algorithm used
% tau    = 0.45;         % BP regularization parameter or OMP sparsity

% method='mli';
% tau = 0.45;
method = 'find_why_low_corr'
tau = [0.2 0.1];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% DATA: RANDOMLY SEPARATED SENSORS, NO MULTIPATH (max M = 240)
%   Note: M = 20 gives good results (prediction correlation coefficient 
%   of greater than 0.99)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% load simdata1_ideal;     % Measured data (same as ideal data)
% load simdata1;           % Ideal data      
% load simdata2_predict;   % Data to predict (same as ideal data right now)
% method = 'omp';          % Optimization algorithm used 
% tau    = 3;              % BP regularization parameter or OMP sparsity
% UNCOMMENT THE NEXT TWO LINES TO USE BASIS PURSUIT DENOISING
% method = 'bp';           % Optimization algorithm used 
% tau    = 0.00001;        % BP regularization parameter or OMP sparsity

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% DATA: LINEARLY ALIGNED SENSORS, DENSE SPACING (0.79 mm) (max M = 250)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% load simdata1_line_ideal;  % Ideal/measure data (arranged in a line)
% load simdata1_line;        % Ideal/measure data (arranged in a line)
% load simdata2_predict;     % Data to predict (placed around the medium)
% tau = 0.00001;             % Basis pursuit regularization parameter
% method = 'omp';            % Optimization algorithm used 
% tau    = 3;                % BP regularization parameter or OMP sparsity
% UNCOMMENT THE NEXT TWO LINES TO USE BASIS PURSUIT DENOISING
% method = 'bp';             % Optimization algorithm used 
% tau    = 0.00001;          % BP regularization parameter or OMP sparsity

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% DATA: LINEARLY ALIGNED SENSORS, SPARSE SPACING (1.6 mm) (max M = 250)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% load simdata2_line_ideal;  % Ideal/measure data (arranged in a line)
% load simdata2_line;        % Ideal/measure data (arranged in a line)
% load simdata2_predict;     % Data to predict (placed around the medium)
% method = 'omp';            % Optimization algorithm used 
% tau    = 3;                % BP regularization parameter or OMP sparsity
% UNCOMMENT THE NEXT TWO LINES TO USE BASIS PURSUIT DENOISING
% method = 'bp';             % Optimization algorithm used 
% tau    = 0.00001;          % BP regularization parameter or OMP sparsity


% NUMBER OF MEASURMENTS USED
M  = 200;

% SET WAVENUMBER DOMAIN PARAMETERS
VN = 500;               % Number of wavenumbers [1/m]
K1 = 2000;              % Maximum wavenumber [1/m]
K0 = K1/VN;             % Minimum wavenumber [1/m]

% SET FREQUENCY DOMAIN PARAMETERS
dF = 1e3;               % Frequency step size [Hz]
F1 = 1000e3;             % Maximum frequency [Hz]
F0 = 10e3;             % Minimum frequency [Hz]

% SPARSE WAVENUMBER ANALYSIS OPTIONS
optMatrix    = false;    % Solve the basis pursuit for all frequencies at once
optNomalized = false;   % Use distortionless constraint
iterPlot     = false;    % Plot result at each frequency


% ******************************************************************* %
% NO NEED TO CHANGE ANYTHING BEYOND THIS POINT
% ******************************************************************* %

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% GET ALL SENSORS AND TRANSMITTER USED 
%   (FOR PLOTTING LATER)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Rx = Rx(:); Tx = Tx(:);
RxSize = cell2mat(cellfun(@size, Rx, 'UniformOutput', false));
for ii = 2:size(RxSize,1), RxSize(ii,:) = RxSize(ii,:) + RxSize(ii-1,:); end
TM = sum(RxSize(:,1) < M)+1;
Rx0 = cell2mat(Rx); Rx0 = Rx0(1:M,:);
Tx0 = cell2mat(Tx(1:TM));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ROUND FREQUENCY INFOMRATION TO FIT DATA
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Q = size(cell2mat(xm),1);
r = floor(Q/2)+1; f = ifftshift((((1:Q)-r)/Q))*Fs; 
[F1x F1n] = min(abs(f-F1)); [F0x F0n] = min(abs(f-F0));
dFn = round(dF/Fs*Q);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SET DOMAINS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fn = (F0n:dFn:F1n).';          % Frequency index space
k = linspace(K0, K1, VN).';    % Sampled wavenumber space
t = 1/Fs:1/Fs:Q/Fs;            % Time space

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% COMPUTE DISTANCES BETWEEN SENSORS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d  = cell2mat(arrayfun(@(ii) dist(Rx{ii} , Tx{ii}.' ), 1:length(Rx) , 'UniformOutput', false).'); 
dp = cell2mat(arrayfun(@(ii) dist(Rxp{ii}, Txp{ii}.'), 1:length(Rxp), 'UniformOutput', false).'); 
dm = cell2mat(arrayfun(@(ii) dist(Rxm{ii}, Txm{ii}.'), 1:length(Rx) , 'UniformOutput', false).'); 
d  = d(1:M);    
dm = dm(1:M);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% COMPUTE FOURIER TRANSFORM OF DATA
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X  = fft(cell2mat(x));         % Fourier transform of ideal data (no corruption)
Xp = fft(cell2mat(xp));        % Fourier transform of data to predict
Xm = fft(cell2mat(xm));        % Fourier transform of measured data 
X  = X(:,1:M);
Xm = Xm(:,1:M);


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PERFORM SPARSE WAVENUMBER ANALYSIS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
V = sparse(VN, Q);
V(:,fn) =  swa_mli_v2_ori( k, dm, Xm(fn,:), tau, 'optMatrix', optMatrix, 'optNomalized', optNomalized, 'plot', iterPlot, 'method', method);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PERFORM SPARSE WAVENUMBER DENOISING
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y = 2*real(sws_mli_v1( k, dm, V ));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PERFORM SPARSE WAVENUMBER SYNTHESIS (PREDICTION)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z = 2*real(sws_mli_v1( k, dp, V ));

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PROCESS MEASURED DATA FOR COMPARISON
%   These lines apply a window to the frequency domain of each 
%   measurement so that we observe the same frequencies as used by 
%   sparse wavenumber analysis
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Xc  = zeros(size(X));
Xmc = zeros(size(Xm));
Xpc = zeros(size(Xp));
Xc(F0n:dFn:F1n,:)  = X(F0n:dFn:F1n,:);  xc = 2*real(ifft(Xc));
Xmc(F0n:dFn:F1n,:) = Xm(F0n:dFn:F1n,:); xmc = 2*real(ifft(Xmc));
Xpc(F0n:dFn:F1n,:) = Xp(F0n:dFn:F1n,:); xpc = 2*real(ifft(Xpc));


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% COMPUTE CORRELATION COEFFICIENTS BETWEEN SIGNALS
% GENERATED BY SWS AND TRUE SIGNALS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ccd = arrayfun(@(i) y(:,i)'*xc(:,i)/norm(y(:,i))/norm(xc(:,i)), 1:size(y,2));
ccs = arrayfun(@(i) z(:,i)'*xpc(:,i)/norm(z(:,i))/norm(xpc(:,i)), 1:size(z,2));

% save('ccd_bp', 'ccd'); 
% save('ccs_bp', 'ccs');
%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% WRITE RESULTS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fprintf('-----------------------------------------------------\n')
fprintf('Denoising results \n')
fprintf('  Average correlation coefficient (/w true result): %f \n', mean(ccd))
fprintf('Prediction results \n')
fprintf('  Average correlation coefficient (/w true result): %f \n', mean(ccs))
fprintf('-----------------------------------------------------\n\n')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PLOT RESULTS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% PLOT SENSOR LOCATIONS
figure(2)
plot(Rx0(:,1), Rx0(:,2), 'LineStyle', 'none', 'MarkerSize', 10, 'Marker', 's', 'color', [37   64  97]/255);
hold on; plot( Tx0(:,1), Tx0(:,2), 'd', 'linewidth', 4, 'MarkerSize', 4, 'color', [99   37  35]/255); hold off;
axis([0 1.22 0 1.22])
legend('Sensors', 'Transmitters','Location','Best')
xlabel('Plate width [m]')
ylabel('Plate length [m]')

% PLOT FREQUENCY-WAVENUMBER REPRESENTATION
figure(3)
imagesc( ifftshift(f)/1000, k, 10*log10(ifftshift(abs(V)/max(max(abs(V))),2)), [-10 0] )
axis([F0/1000 F1/1000 K0 K1]); axis xy;
xlabel('Frequency [kHz]');
ylabel('Wavenumber [m^{-1}]');

% PLOT SPARSE WAVENUMBER DENOISING RESULTS
figure(4)
subplot(311); plot(t*1000, xmc(:,2));
text(max(get(gca, 'xlim'))*.975, max(get(gca, 'ylim')*.875), 'Measured signal with multipath', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
xlabel('Time [ms]'); ylabel('Amplitude');
subplot(312); plot(t*1000, y(:,2));
text(max(get(gca, 'xlim'))*.975, max(get(gca, 'ylim')*.875), 'Sparse wavenumber denoised signal', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
xlabel('Time [ms]'); ylabel('Amplitude')
subplot(313); plot(t*1000, xc(:,2));
text(max(get(gca, 'xlim'))*.975, max(get(gca, 'ylim')*.875), 'True signal without multipath', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
xlabel('Time [ms]'); ylabel('Amplitude')

% PLOT SPARSE WAVENUMBER SYNTHESIS RESULTS
figure(5)
subplot(211); plot(t*1000, z(:,2)); 
text(max(get(gca, 'xlim'))*.975, max(get(gca, 'ylim')*.875), 'Sparse wavenumber synthesis predicted signal', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
xlabel('Time [ms]'); ylabel('Amplitude');
subplot(212); plot(t*1000, xpc(:,2));
text(max(get(gca, 'xlim'))*.975, max(get(gca, 'ylim')*.875), 'True signal', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
xlabel('Time [ms]'); ylabel('Amplitude');

% PLOT CORRELATION COEFFICIENTS
figure(6)
subplot(211); plot([1:length(ccd)].', ccd, 'linewidth', 2); 
title('Correlation coefficients between true and denoised signals')
xlabel('Measurement Number [no order]'); ylabel('Correlation Coefficient');
subplot(212); plot([1:length(ccs)].', ccs, 'linewidth', 2);
title('Correlation coefficients between true and predicted signals')
xlabel('Measurement Number [no order]'); ylabel('Correlation Coefficient');

