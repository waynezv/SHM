function [IDX] = AsymSpectral(CKSym,n,mode)

N = size(CKSym,1);
MAXiter = 1000; % Maximum iteration for KMeans Algorithm
REPlic = 100; % Replication for KMeans Algorithm

if mode == 1
    % Method 1: Unnormalized Method
    DKU = diag( sum(CKSym) );
    LapKU = DKU - CKSym;
    [uKU,sKU,vKU] = svd(LapKU);
    f = size(vKU,2);
    kerKU = vKU(:,f-n+1:f);
    svalKU = diag(sKU);
    IDX = kmeans(kerKU,n,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
    
elseif mode==2
    % Method 2: Random Walk Method
    DKN=( diag( sum(CKSym) ) )^(-1);
    LapKN = speye(N) - DKN * CKSym;
    [uKN,sKN,vKN] = svd(LapKN);
    f = size(vKN,2);
    kerKN = vKN(:,f-n+1:f);
    svalKN = diag(sKN);
    IDX = kmeans(kerKN,n,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
    
elseif mode==3
    % Method 3: Normalized Symmetric
    DKS = ( diag( sum(CKSym) ) )^(-1/2);
    LapKS = speye(N) - DKS * CKSym * DKS;
    [uKS,sKS,vKS] = svd(LapKS);
    f = size(vKS,2);
    kerKS = vKS(:,f-n+1:f);
    for i = 1:N
        kerKS(i,:) = kerKS(i,:) ./ norm(kerKS(i,:));
    end
    svalKS = diag(sKS);
    IDX = kmeans(kerKS,n,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
    
elseif mode==4
    % Method 4: Total variation - Symmetric
    LapTVS = (speye(N) - CKSym)'*(speye(N) - CKSym);
    [uTVS,sTVS,vTVS] = svd(LapTVS);
    f = size(vTVS,2);
    kerTVS = vTVS(:,f-n+1:f);
    for i = 1:N
        kerTVS(i,:) = kerTVS(i,:) ./ norm(kerTVS(i,:));
    end
    svalTVS = diag(sTVS);
    IDX = kmeans(kerTVS,n,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
    
elseif mode==5
    % Method 5: Total variation - Asymmetric
    P = diag(sum( CKSym, 2).^(-1))*CKSym;
    LapTVA = (speye(N) - P)'*(speye(N) - P);
    [uTVA,sTVA,vTVA] = svd(LapTVA);
    f = size(vTVA,2);
    kerTVA = vTVA(:,f-n+1:f);
    for i = 1:N
        kerTVA(i,:) = kerTVA(i,:) ./ norm(kerTVA(i,:));
    end
    svalTVA = diag(sTVA);
    IDX = kmeans(kerTVA,n,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
else
    display('Error! Please indicate mode between 1 and 5')
    return
end

%SingVals = [svalKU,svalKN,svalKS,svalTVS,svalTVA];
%LapKernel(:,:,1) = kerKU;
%LapKernel(:,:,2) = kerKN;
%LapKernel(:,:,3) = kerKS;
%LapKernel(:,:,4) = kerTVS;
%LapKernel(:,:,5) = kerTVA;

