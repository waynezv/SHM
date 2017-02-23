Vn = magic(7)

orimap=abs(Vn)>0.02;
orimap

            feat=[];
%               feat = zeros(size(orimap));
            for i=1:size(orimap,1)
                for j=1:size(orimap,2)
                    if(orimap(i,j)>0)
                        feat=[feat;[i,j]];
                    end
                end
            end
            
            feat
            
            feat2=[feat(:,1) feat(:,2)-feat(:,1)];
            
              x=feat2';
               [a, model, L] = vbgm(x, 6);
               
               label_name=unique(a);
               kmean_K=length(label_name);
               
              xx = bsxfun(@times, magic(5),eye(5));
              xlswrite('test_xx.xls', xx);
              yy = magic(25);
              xlswrite('test_yy.xls', yy);
              
              i = 1;
               fprintf('ss');
               fprintf('%d',i);
               fprintf('ss\n');
               
              
              fprintf('%s',i,'to xls...',\n);
              
              disp(['varargin{' num2str(v) '} class is ' class(varargin{v})]);
            
            
            
            