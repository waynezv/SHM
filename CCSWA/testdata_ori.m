% ---------------------------------------------
% OPTIONS
% ---------------------------------------------
fn = 10:1000;

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

%V = swa( k, d, X(200:800,:), .4, 'method', 'bp', 'plot', true );
V = swa( k, d, X(fn,:), 2, 'method', 'omp' );
save VrealDataOri V;
%sparse wavenumber synthesis
z = 2*real(sws_mli_v1( k, d, V ));
% save swsOri z;
%%
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