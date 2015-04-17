% ---------------------------------------------
% OPTIONS
% ---------------------------------------------
fn = 10:100;

% ---------------------------------------------
% GET DATA
% ---------------------------------------------
[meta, cfg] = extract_data(@data_20121029_localize_baseline, 1:34);


%%
% ---------------------------------------------
% GET DISPERSION CURVES
% ---------------------------------------------

k = (2:2:1000).';
Tx = arrayfun(@(ii) (meta.Tx{ii}.'*ones(1,size(meta.Rx{ii},1))).', 1:numel(meta.Rx), 'UniformOutput', false);
%d = diag(dist2(cell2mat(meta.Rx), cell2mat(Tx.').'));
d = (dist2(cell2mat(meta.Rx), cell2mat(Tx.').'));

x = cell2mat(meta.x.');
x(1:110,:) = 0; % Remove cross-talk
x = velwindow(x, d, 2000/meta.Fs);
X = fft(x);
% 此处不知Q如何设置
Q = size(x,1);

%V = swa( k, d, X(200:800,:), .4, 'method', 'bp', 'plot', true );
% V = swa( k, d, X(fn,:), 2, 'method', 'omp' );
%%
method='mli';
tau    = 0.45;

% SPARSE WAVENUMBER ANALYSIS OPTIONS
optMatrix    = false;    % Solve the basis pursuit for all frequencies at once
optNomalized = false;   % Use distortionless constraint
iterPlot     = false;    % Plot result at each frequency


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PERFORM SPARSE WAVENUMBER ANALYSIS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%此处sparse error: Sparse matrix sizes must be non-negative integers less than MAXSIZE as defined by COMPUTER.
% V = sparse(V, Q);


    V(:,fn) =  swa_mli_md_1( k, d, X(fn,:), tau, 'optMatrix', optMatrix, 'optNomalized', optNomalized, 'plot', iterPlot, 'method', method);
save V_real V;
%%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PERFORM SPARSE WAVENUMBER DENOISING
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y = 2*real(sws_mli_v1( k, d, V ));
% ---------------------------------------------
% PLOT
% ---------------------------------------------

figure(1)
N = meta.N; Fs = meta.Fs; t = 1/Fs:1/Fs:N/Fs;
plot(t, x)
xlabel('Time')
ylabel('Voltage')


figure(2)
f = linspace(0, Fs, N);
imagesc(f(fn)/1000, k, 10*log10(abs(V)/max(max(abs(V)))), [-25 0])
axis xy;
ylabel('Wavenumber')
xlabel('Frequency [kHz]')
colormap(fireprint)
a = colorbar;
xlabel(a, '[dB]')