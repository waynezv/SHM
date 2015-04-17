clear all;
load('example.mat');
V(:,fn) =  swa_mli_v1( k, dm, Xm(fn,:), tau, 'optMatrix', optMatrix, 'optNomalized', optNomalized, 'plot', iterPlot, 'method', method);
figure(36)
imagesc( ifftshift(f)/1000, k, 10*log10(ifftshift(abs(V)/max(max(abs(V))),2)), [-10 0] )
axis([F0/1000 F1/1000 K0 K1]); axis xy;
xlabel('Frequency [kHz]');
ylabel('Wavenumber [m^{-1}]');