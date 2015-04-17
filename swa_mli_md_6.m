function [V Phi] = swa_mli_md_6( k, d, X, tau, varargin )
% todo：
%SWA  Sparse wavenumber analysis
%   [V Phi DrADk] = SWA( K, D, X, TAU, OPTIONS ) transforms wave data into
%   a sparse wavenumber domain. 
%
%   INPUTS:   
%       k: An N-by-1 vector of wavenumbers
%       d: An M-by-1 vector of distancess
%       X: An Q-by-M cell or matrix of data
%     tau: Regularization parameter for optimization basis pursuit 
%          denoising or sparsity value for orthogonal matching pursuit
%  
%   OPTIONS: 
%           'method': 'bp' for basis pursuit denoising, 'omp' for
%                     orthogonal matching pursuit
%        'optMatrix': Solve the basis pursuit for all frequencies at once
%    'optNormalized': If true, use a distortionless constrant
%             'plot': If true, plot result after each frequency
%
%   OUTPUTS:
%     V: Sparse wavenumber solution
%     Phi: Wave propagation frame used
%
%   see also: sws
%
% -------------------------------------------------------------------------
% Code: Written by: Joel B. Harley
% Modified by: Ming Li
% last updated: July 13th 2014
% -------------------------------------------------------------------------
%
    
    % --------------------------------------------------------------------
    % MANAGE INPUT ARGUMENTS
    % --------------------------------------------------------------------

    % CHECK NUMBER OF ARGUMENTS
    if nargin < 4, error('SWA requires 4 or more input arguments.'); end 

    % FIX ARGUMENT FORMATS
    if  iscell(X), X = cell2mat(X); end     % Make a matrix
    %d = d(:);                              % Make a column vector
    %k = k(:);                              % Make a column vector
    
    
    % SET DEFUALT OPTIONAL ARGUMENTS
    opt.optMatrix     = false;              % Use fully matrix optimization (memory intensive) 
    opt.optNomalized  = false;              % Use normalized optimization
    opt.plot          = false;              % Plot results as they run (doesn't work if optMatrix = true)
    opt.method        = 'mli';               % Optimization method
    
    % PARSE ARGUMENTS AND SIMPLIFY SOME ARGUMENT NAMES
    if ~isempty(varargin), opt = parseArgs(opt, varargin{:}); end
    
    % --------------------------------------------------------------------
    % PERFORM FUNCTION
    % --------------------------------------------------------------------
    penalty1=1;
    penalty2=0.5;
    % DEFINE LENGTHS
    Q = size(X, 1);                         % Number of frequencies
%     Q
    % NORMALIZE INPUT DATA
    EX = sqrt(sum(abs(X).^2,2));            % Data normalization factors (at each frequency)
    Xn = bsxfun(@times, X, 1./EX);          % Normalized data
%     about Xn
    % BUILD PROPOGATION MATRICES
    [Phi Dk rho] = waveframe(d, k);                  % Compute propogation frame
%     about Phi 

    % PERFORM OPTIMIZATION TO COMPUTE SPARSE WAVENUMBER SOLUTION
    switch opt.method
        case 'bp'
            if   opt.optMatrix, Vn = bp(Phi, Xn.', tau, opt.optNomalized);  
            else Vn = cell2mat(arrayfun(@(ii) countLoop(@bp, ii, Q, opt.plot*k, Phi, permute(Xn(ii,:,:), [2 1 3]), tau, opt.optNomalized), 1:Q, 'UniformOutput', false )); end
        case 'omp'
            if   opt.optMatrix, Vn = omp(Phi, Xn.', tau, opt.optNomalized);  
            else Vn = cell2mat(arrayfun(@(ii) countLoop(@omp, ii, Q, opt.plot*k, Phi, permute(Xn(ii,:,:), [2 1 3]), tau, opt.optNomalized), 1:Q, 'UniformOutput', false )); end
        case 'mli'       
            weight=ones(size(Phi,2),size(Xn,1));               
       load Vn_unweighted_md3;
                
                orimap=abs(Vn)>0.02;% 此处取了绝对值 可否用complex？
                feat = [];
                for i=1:size(orimap,1)
                    for j=1:size(orimap,2)
                        if(orimap(i,j)>0)
                            feat=[feat;[i,j]];
                        end
                    end
                end
               % [U,S,V] = svd(feat);
               % feat2=feat*V;
               % 此处存疑
               feat2=[feat(:,1) feat(:,2)-feat(:,1)];
    %            kmean_K=3;
    %             [a b]=kmeans(feat2,kmean_K);
    
               x=feat2';
               [a, ~] = vbgm(x, 10);
               label_name=unique(a);
               kmean_K=length(label_name);% 聚类数
               % 拟合
               for k=1:kmean_K
                   % 二项式拟合系数
                [P(k,:)]= polyfit(feat(find(a==k),2),feat(find(a==k),1),2);% 其他插值方法？
%                 [P(k,:)]= interp1(feat(find(a==k),2),feat(find(a==k),1),2);
                % 二项式拟合在feat各点处的取值
                t1=polyval(P(k,:),feat(find(a==k),2));
                errorpolyfit(k)=mean((feat(find(a==k),1)-t1).^2);% Square Mean Error
               end               
% 以下计算新的weight 来找回缺失的outlier
               for k=1: kmean_K
                   % 如果和拟合曲线误差小于50 找回这些点
                       if(errorpolyfit(k)<50)
                           for i=1:Q
                               if(i>min(feat(find(a==k),2)))
                                  t1 = P(k,1)*i^2 + P(k,2)*i +P(k,3);
                                  t2 = round(t1);
                                  %t3=max(t2-1,1);
                                  %t4=min(t2+1,size(weight,1));
                                  %weight(t3:t4,i)=weight(t3:t4,i)*penalty2; 
                                  % 此处存疑
                                    if(((t2-1)>=1)&&((t2+1)<size(weight,1)))
                                        weight(t2-1:t2+1,i)=penalty2;
                                    end
                              end
                           end
                       end
               end
                fprintf('saving weight after ploy...\n');
                save weight_fit_md5 weight;
    %             matlabpool open 12;
                   parfor j=1:Q
                         Vn(:,j) = bp2(Phi, Xn(j,:).', weight(:,j),tau, opt.optNomalized);  
                         j
                   end
                   fprintf('saving Vn after poly...\n');
                   save Vn_poly_real_md5 Vn;
    %                matlabpool close;         
%             cvx_end
            [BW,thresh,gv,gh] = edge(abs(Vn),'sobel',0.02);% sobel 改进，0.02？
                        for j=2:Q-1
                            for k=2:(size(Phi,2)-1)
                             if((abs(Vn(k,j))>0.001)&&(sum(sum(BW(k-1:k+1,j-1:j+1)))==0))
                                %weight(k,j)=weight(k,j)*penalty1;
                                weight(k,j)=penalty1;% 改进weight？variational
                             end
                            end
                        end
%                 fprintf('saving weight...\n');
%                 save weight_real weight;
               
    end
    
    V=rho.*inv(Dk)*Vn*diag(EX);
    % DEBIAS RESULTS
    Y   = ((Phi*Dk*V)./rho).';                      % Generate denoised signal
    mu  = cell2mat(arrayfun(@(ii) X(ii,:)*Y(ii,:)' ./ norm(Y(ii,:))^2, 1:Q, 'UniformOutput', false ));
    V = bsxfun(@times, mu, V);          % Frequency-wavenumber representation
    
end


function [Phi Dk rho] = waveframe(d, k) 
%WAVEFRAME  Compute a wave propagation frame
%   [Phi Dk] = WAVEFRAME(D, A, K)  computes the a matrix describing wave
%   propagation for a collection of chosen wavenumbers and distances
%
%   INPUTS: 
%       d: An M-by-1 vector of distancess
%       a: An M-by-1 vector of amplitudes 
%       k: An N-by-1 vector of wavenumbers 
%
%  OUTPUTS: 
%     Phi: An M-by-N basis matrix. Each column represents an impulse 
%          response for a given wavenumber and set of distances
%     rho: Normalization energy
%      Dk: An N-by-N diagonal matrix of wavenumber weights
%
% -------------------------------------------------------------------------
% Code: Written by: Joel B. Harley
% Last updated: May 5, 2014
% -------------------------------------------------------------------------
%

    % INITIALIZE LENGTHS
    M = size(d,1); N = length(k);

    % COMPUTE THE PROPAGATION MATIX PHI
    A = exp(-1j*d*k.');                 % Wave propagation frame
    Dr = spdiags(1./sqrt(d), 0, M, M);  % Distance attenuation
    Dk = spdiags(1./sqrt(k), 0, N, N);  % Wavenumber attenutation
    rho = 1./norm(Dr,'fro');            % Normalization: 1./sqrt(sum(1./d))
    Phi = rho*Dr*A;                     % Progagation frame
end


function V = bp( A, X, tau, normFlg ) 
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
        
        
        
       minimize( cst + tau*sum(norms( V , 1, 1 )) );
        if normFlg == true
            subject to   % Distortionless constraint
            for p = 1:P, trace(X(:,:,p)'*(A(:,:,p)*V))/norm(X(:,:,p), 'fro') == 1;  end
        end    
    cvx_end


end



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
        minimize( cst + tau*sum(norms( V.*weight, 1, 1 )) );
        if normFlg == true
            subject to   % Distortionless constraint
            for p = 1:P, trace(X(:,:,p)'*(A(:,:,p)*V))/norm(X(:,:,p), 'fro') == 1;  end
        end    
    cvx_end


end


function x = omp( A, b, K, normFlg ) 
%myomp  Compute the sparse basis pursuit solution
%   Z = BP(PSI, X, TAU, NORMFLAG)  computes the sparse solution
%   Z given a basis matrix PSI and observations X.
%
%   INPUTS: 
%       A: An M-by-N basis matrix
%       b: An M-by-1 vector of observations
%       K: Number of sparse components
% normFlg: Not used -- included to match form with BP function
%
%  OUTPUTS: 
%       x: The N-by-1 wavenumber basis pursuit denoising result
%
% -------------------------------------------------------------------------
% Code: Written by: Joel B. Harley
% Last updated: May 5, 2015 
% -------------------------------------------------------------------------
%
    
    % DEFINE SIZES 
    N   = size(A,2);            % Number of atoms
    M   = size(A,1);            % Aize of atoms
    Q   = size(b,2);            % Number of outputs/frequencies

    % INITIALIZE VARIABLES
    x         = zeros(N,Q);     % Solution
    x_T       = zeros(K,Q);     % Solution with only non-zero indices
    indx_set  = zeros(K,Q);     % Indices with non-zero values
    atoms     = zeros(M,K,Q);   % Chosen dictionary atoms for each frequency

    % INIATIVE ALGORITHM 
    r   = b;                    % Initial residual
    xr  = A'*r;                 % Initial solution from residual xr
    
    % LOOP OVER NUMBER OF SPARSE COMPONENTS
    for k = 1:K
        
        % ----------------------------------------------------------- %
        % -- FIND CORRESPONDING INDICES
        % ----------------------------------------------------------- %
        [~,ind_new]   = max(abs(xr), [], 1);
        indx_set(k,:) = ind_new;
        atoms(:,k,:)  = permute(A(:,ind_new), [1 3 2]);

        % ----------------------------------------------------------- %
        % -- UPDATE RESIDUAL
        % ----------------------------------------------------------- %
        % LOOP OVER OUTPUTS
        for q = 1:Q
            x_T(1:k,q) = atoms(:,1:k,q) \ b(:,q);         % Find least-squares fit
            x( indx_set(1:k,q), q )   = x_T(1:k,q);       % Places results in full vector
            r(:,q)  = b(:,q) - atoms(:,1:k,q)*x_T(1:k,q); % Find new residual
        end
        
        % ----------------------------------------------------------- %
        % -- COMPUTE SOLUTION FROM RESIDUAL xr FOR NEXT ITERATION
        % ----------------------------------------------------------- %
        if k < K, xr  = A'*r; end

    end

end

function varargout = countLoop( fnc, n, N, k, varargin )
%COUNTLOOP  Displays the time till complete of a function loop 
%   VARARGOUT = COUNTLOOP( FNC, n, N, k, VARARGIN ) runs the function 
%   FNC, displays the loop count n/N, and displays the estimated time to
%   completion
%
%   INPUTS:   
%       FNC: A function handle to run
%         n: Current interation in loop
%         N: Last interation in loop
%         K: If non-zero, plots the function output with "K" defined as 
%            the horizontal axis
%  VARARGIN: Input parameters to the function "fnc"
%
%   OUTPUTS:
% VARARGOUT: Output of function "fnc"
%
% -------------------------------------------------------------------------
% Code: Written by: Joel Harley
% Last updated: November 7, 2013 
% -------------------------------------------------------------------------
%

    % DISPLAY TIME INFORMATION
    global tm; if isempty(tm), tm = 0; end
    fprintf('%08i / %08i [Time left: %s]', n, N, datestr(tm/24/3600*(N-n+1), 'HH:MM:SS')); ts = tic; 
    
    % RUN FUNCTION
    M = nargout(fnc); if M == -1, M = 1; end;
    [varargout{1:M}] = fnc(varargin{:});

    % REFRESH TIME INFORMATION
    fprintf(repmat('\b', 1, 41)); 
    
    % REFRESH TIME INFORMATION
    tm = (toc(ts) + tm)/2;
    
    % PLOT RESULT IF REQUESTED
    if any(k)
        if length(varargout{1}) == length(k)
            figure(1);
            plot(k, abs(varargout{1})); axis xy;
            xlabel('Wavenumber [m^{-1}]')
            drawnow;
        elseif length(varargout{1}) == (length(k)^2)
            figure(1)
            V0 = reshape(varargout{1}, length(k), length(k));
            imagesc(k,k,abs(squeeze(V0(:, :))))
            colormap(copper)
            drawnow;
        else
            figure(1);
            plot(abs(varargout{1})); axis xy;
            xlabel('Wavenumber [samples]')
            drawnow;
        end
    end
    
end



function options = parseArgs(options, varargin)

    % ---------------------------------------------------
    % CODE TAKEN FROM STACKOVERFLOW
    %   http://stackoverflow.com/questions/2775263/how-to-deal-with-name-value-pairs-of-function-arguments-in-matlab
    % ---------------------------------------------------
    
    % GET OPTIONS NAMES
    optionNames = fieldnames(options);

    % COUNT ARGUMENTS
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('SWA needs propertyName/propertyValue pairs')
    end
    
    % PARSE ARGUMENTS
    for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
       inpName = pair{1}; % make case insensitive

       if any(strcmpi(inpName,optionNames))
          % overwrite options. If you want you can test for the right class here
          % Also, if you find out that there is an option you keep getting wrong,
          % you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
          options.(inpName) = pair{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end

end

