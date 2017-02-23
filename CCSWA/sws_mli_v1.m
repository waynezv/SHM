function x = sws_mli_v1( k, d, V )
%SWA  Sparse wavenumber analysis
%   [VS Y A Dr] = SWA( RX, TX, RTA, X, OPTIONS ) transforms wave data into
%   a sparse wavenumber domain. 
%
%   INPUTS:       
%       k: An N-by-1 vector of wavenumbers
%       d: An M-by-1 vector of distancess
%       V: Sparse wavenumber solution from SWA
%
%   OUTPUTS:
%       x: Synthesized data
%
% -------------------------------------------------------------------------
% Code: Written by: Joel B. Harley
% Last updated: November 7, 2013 
% -------------------------------------------------------------------------
%

    % CHECK INPUTS 
    error(nargchk(3, 3, nargin));
    
    % FORCE COLUMN VECTORS
    k = k(:);
    
    % DEFINE LENGTHS
    M = size(d, 1);           % Number of signals to synthesize
    N = size(k, 1);           % Number of wavenumbers
    
    
    if size(d,2) == 1                  % COMPUTE PROPAGATION FRAME
        A = exp(-1j*d*k.');                % Wave propogation frame
        Dr = spdiags(1./sqrt(d), 0, M, M); % Distance attenuation
        Dk = spdiags(1./sqrt(k), 0, N, N); % Wavenumber attenutation 
                                           % Note: Dk is assumed to be part of
                                           %       V and so not used here
        Phi = 1./norm(Dr, 'fro')*Dr*A;
        DrADk=Dr*A*Dk;
      else
        error('waveframe:sws:invalidDimensions', 'd can only be M-by-1.') 
    end
    
        
    % SYNTHESIZE SIGNALS
    x = ifft((DrADk*V).');  
        
end

