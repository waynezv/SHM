function V = bp2( A, X, weight,tau, normFlg ) 
%BP  Compute the sparse basis pursuit solution
%   Z = BP(PSI, X, TAU, NORMFLAG)  computes the sparse solution
%   Z given a basis matrix PSI and observations X.
%
%   INPUTS: 
%       A: An M-by-N basis matrix
%       X: An M-by-1 vector of observations
%     tau: Regularization parameter for optimization
% normFlg: If true, use a distortionless constrant
%
%  OUTPUTS: 
%       V: The N-by-1 wavenumber basis pursuit denoising result
%
% -------------------------------------------------------------------------
% Code: Written by: Joel B. Harley
% Last updated: November 7, 2013 
% -------------------------------------------------------------------------
%
 
    % INITIALIZE VARIABLES
    Q = size(X,2);      % Number of frequencies
    N = size(A,2);      % Number of wavenumbers
    P = size(A,3);      % Number of pathes
    cst = 0;            % Initialize cost function

    % RUN CVX OPTIMIZATION
    cvx_begin quiet
        variable V(N,Q) complex

        % CONSTRUCT COST FUNCTION
        for p = 1:P, cst = cst + pow_pos(norm( A(:,:,p)*V - X(:,:,p) , 'fro' ),2); end

        % MINIMIZE COST 
        
        
        minimize( cst + tau*sum(norms( V.*weight , 1, 1 )) );
        if normFlg == true
            subject to   % Distortionless constraint
            for p = 1:P, trace(X(:,:,p)'*(A(:,:,p)*V))/norm(X(:,:,p), 'fro') == 1;  end
        end    
    cvx_end
end
