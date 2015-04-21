function [V Phi] = swa_mli_v2( k, d, X, tau, varargin )
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
    % DEFINE LENGTHS
    Q = size(X, 1);                         % Number of frequencies
    
    % NORMALIZE INPUT DATA
    EX = sqrt(sum(abs(X).^2,2));            % Data normalization factors (at each frequency)
    Xn = bsxfun(@times, X, 1./EX);          % Normalized data
    
    % BUILD PROPOGATION MATRICES
    [Phi Dk rho] = waveframe(d, k);                  % Compute propogation frame
    
    % PERFORM OPTIMIZATION TO COMPUTE SPARSE WAVENUMBER SOLUTION
    switch opt.method
        case 'bp'
            if   opt.optMatrix, Vn = bp(Phi, Xn.', tau, opt.optNomalized);  
            else Vn = cell2mat(arrayfun(@(ii) countLoop(@bp, ii, Q, opt.plot*k, Phi, permute(Xn(ii,:,:), [2 1 3]), tau, opt.optNomalized), 1:Q, 'UniformOutput', false )); end
        case 'omp'
            if   opt.optMatrix, Vn = omp(Phi, Xn.', tau, opt.optNomalized);  
            else Vn = cell2mat(arrayfun(@(ii) countLoop(@omp, ii, Q, opt.plot*k, Phi, permute(Xn(ii,:,:), [2 1 3]), tau, opt.optNomalized), 1:Q, 'UniformOutput', false )); end
        case 'mli'
<<<<<<< HEAD
            % LOAD DATA
            addpath ./data_mat
            load Vn_1m_li.mat;
            % INITIALIZE
            nClus = 35; % set # of clusters
            penalty1=1.5;
            penalty2=0.5;
            
%              % INITIALIZE WEIGHT
%              weight=ones(size(Phi,2),size(Xn,1));
%             for i=1:2
%                 if i==1
%                      weight=ones(size(Phi,2),size(Xn,1));
%                 else
%                     [BW,thresh,gv,gh] = edge(abs(Vn),'sobel',0.02);% sobel ????????????0.02????
%                         for j=2:Q-1
%                             for k=2:(size(Phi,2)-1)
%                              if((abs(Vn(k,j))>0.001)&&(sum(sum(BW(k-1:k+1,j-1:j+1)))==0))
%                                 %weight(k,j)=weight(k,j)*penalty1;
%                                 weight(k,j)=penalty1;% ????????weight????variational
%                              end
%                             end
%                         end
%                 end
%             % COMPUTE VN WITH EDGE WEIGHTING
%             matlabpool open 12
%                 parfor j=1:Q
%                      Vn(:,j) = bp2(Phi, Xn(j,:).', weight(:,j),tau, opt.optNomalized);  
%                      j
%                 end
%                 fprintf('saving Vn after edging...\n');
% %                 save VnReal105 Vn;
%             matlabpool close
         
            %% EXTRACT FEATURES
=======
             % INITIALIZE WEIGHT
             weight=ones(size(Phi,2),size(Xn,1));
            for i=1:2
                if i==1
                     weight=ones(size(Phi,2),size(Xn,1));
                else
                    [BW,thresh,gv,gh] = edge(abs(Vn),'sobel',0.02);% sobel ????????????0.02????
                        for j=2:Q-1
                            for k=2:(size(Phi,2)-1)
                             if((abs(Vn(k,j))>0.001)&&(sum(sum(BW(k-1:k+1,j-1:j+1)))==0))
                                %weight(k,j)=weight(k,j)*penalty1;
                                weight(k,j)=penalty1;% ????????weight????variational
                             end
                            end
                        end
                end
            % COMPUTE VN WITH EDGE WEIGHTING
            matlabpool open 12
                parfor j=1:Q
                     Vn(:,j) = bp2(Phi, Xn(j,:).', weight(:,j),tau, opt.optNomalized);  
                     j
                end
                fprintf('saving Vn after edging...\n');
                save VnReal105 Vn;
            matlabpool close
            % SUBSTRACT FEATURES
           for ii=1:1
>>>>>>> 069b1039758e285b8de0705197a7ca706af3a154
            orimap=abs(Vn)>0.02;
            feat=[];
            for i=1:size(orimap,1)
                for j=1:size(orimap,2)
                    if(orimap(i,j)>0)
                        feat=[feat;[i,j]];% DIM (i * j) * 2
                    end
                end
            end
           % [U,S,V] = svd(feat);
           % feat2=feat*V;
           feat2=[feat(:,1) feat(:,2)-feat(:,1)]; % same DIM with feat
           % kmean_K=3;
           % [a b]=kmeans(feat2,kmean_K);
           x=feat2'; % DIM 2 * ( i * j )
<<<<<<< HEAD
           
           %% CLUSTER
           % ADJUST NO. OF CLUSTERS TO nClus FOR 1MHz CASE
           [a, model, L] = vbgm_wz_1(x, nClus); % DIM ( i* j ) * 1 
           label_name=unique(a); % Get the label name without repetition
           kmean_K=length(label_name);% Get the number of labels  
           
            % REDUCE CLUSTER BY REMOVING CLUSTERS WITH FEW DOTS
            ClusRed = []; % init reduced # of clusters
            for k = 1 : kmean_K
                nDot = length(find(a==k));
                if  nDot > 15 % assume a cluster size threshold
                    ClusRed = [ClusRed k]; % record cluster label that larger than threshold size
                end
            end

            % ADJUST WEIGHT BY POLYNIMIAL FIT
            for k = ClusRed
                xx = feat(find(a==k),2);
                yy = feat(find(a==k),1);
                [p(k,:),~,mu] = polyfit(xx,yy,2);
                for i = 1:Q
                    if i>min(xx) && i<max(xx)
                         valuePred = polyval(p(k,:),i,[],mu);  
                         weight(round(valuePred), i) = penalty2;
                    end
                end
            end
            
%            % POLYNOMIAL FIT
%            for k=1:kmean_K
%             P(k,:) = polyfit(feat(find(a==k),2),feat(find(a==k),1),2);
%             t1=polyval(P(k,:),feat(find(a==k),2));
%             errorpolyfit(k)=mean((feat(find(a==k),1)-t1).^2);% Square Mean Error
%            end
%            % ADJUST WEIGHT WITH FIT 
%             for k=1: kmean_K
%                % ????????????????????????????????????????????50 ????????????????????
%                    if(errorpolyfit(k)<50)
%                        for i=1:Q
%                            if(i>min(feat(find(a==k),2)))
%                               t1 = P(k,1)*i^2 + P(k,2)*i +P(k,3);
%                               t2 = round(t1);
%                               %t3=max(t2-1,1);
%                               %t4=min(t2+1,size(weight,1));
%                               %weight(t3:t4,i)=weight(t3:t4,i)*penalty2; 
%                               % ????????????????
%                                 if(((t2-1)>=1)&&((t2+1)<size(weight,1)))
%                                     weight(t2-1:t2+1,i)=penalty2;
%                                 end
%                           end
%                        end
%                    end
%            end

           %% COMPUTE VN AGAIN WITH CLUSTER WEIGHTING
            matlabpool open
=======
           % CLUSTER
           % ADJUST NO. OF CLUSTERS TO 4 FOR 1MHz CASE
           [a, model, L] = vbgm_wz_1(x, 4); % DIM ( i* j ) * 1 
           label_name=unique(a); % Get the label name without repetition
           kmean_K=length(label_name);% Get the number of labels           
           % POLYNOMIAL FIT
           for k=1:kmean_K
            P(k,:) = polyfit(feat(find(a==k),2),feat(find(a==k),1),2);
            t1=polyval(P(k,:),feat(find(a==k),2));
            errorpolyfit(k)=mean((feat(find(a==k),1)-t1).^2);% Square Mean Error
           end
           % ADJUST WEIGHT WITH FIT 
            for i=1:Q
                  for k=1: kmean_K
                      if(i>min(feat(find(a==k),2)))
                            t1= P(k,1)*i^2 + P(k,2)*i +P(k,3);
                            t2=round(t1);
%                           weight(t2-1:t2+1,i)=weight(t2-1:t2+1,i)*penalty2; 
                            weight(t2:t2+1,i)=weight(t2:t2+1,i)*penalty2; 
                      end
                  end
            end
           % USE BELOW FOR ALTERNATIVE
            for k=1: kmean_K
               % ????????????????????????????????????????????50 ????????????????????
                   if(errorpolyfit(k)<50)
                       for i=1:Q
                           if(i>min(feat(find(a==k),2)))
                              t1 = P(k,1)*i^2 + P(k,2)*i +P(k,3);
                              t2 = round(t1);
                              %t3=max(t2-1,1);
                              %t4=min(t2+1,size(weight,1));
                              %weight(t3:t4,i)=weight(t3:t4,i)*penalty2; 
                              % ????????????????
                                if(((t2-1)>=1)&&((t2+1)<size(weight,1)))
                                    weight(t2-1:t2+1,i)=penalty2;
                                end
                          end
                       end
                   end
           end
           % COMPUTE VN AGAIN WITH CLUSTER WEIGHTING
            matlabpool open 12
>>>>>>> 069b1039758e285b8de0705197a7ca706af3a154
               parfor j=1:Q
                     Vn(:,j) = bp2(Phi, Xn(j,:).', weight(:,j),tau, opt.optNomalized);  
                     j
               end
            matlabpool close
<<<<<<< HEAD
            printf('saving Vn...\n');
            save('./data_mat/Vn0421.mat', 'Vn');
=======
            
           end
>>>>>>> 069b1039758e285b8de0705197a7ca706af3a154
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
        % wayne in: par
        
        parfor p = 1:P, cst = cst + pow_pos(norm( A(:,:,p)*V - X(:,:,p) , 'fro' ),2); end

        % MINIMIZE COST 
        
        
        minimize( cst + tau*sum(norms( V.*weight , 1, 1 )) );
        if normFlg == true
            subject to   % Distortionless constraint
            parfor p = 1:P, trace(X(:,:,p)'*(A(:,:,p)*V))/norm(X(:,:,p), 'fro') == 1;  end
        end    
    cvx_end
    
    % wayne out


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

