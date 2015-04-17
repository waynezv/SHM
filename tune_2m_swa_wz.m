addpath .\121318
load feat
load feat2   
load weight

Q = 1991;
penalty1=1.5;
penalty2=0.5;
opt.optNomalized  = false;
                
x=feat2';
               [a, model, L] = vbgm(x, 10);
               label_name=unique(a);
               kmean_K=length(label_name);

               for k=1:kmean_K
                [P(k,:)]= polyfit(feat(find(a==k),2),feat(find(a==k),1),2);
                t1=polyval(P(k,:),feat(find(a==k),2));
                errorpolyfit(k)=mean((feat(find(a==k),1)-t1).^2);
               end

               %calculate error

               for k=1: kmean_K

                       if(errorpolyfit(k)<50)
                       for i=1:Q

                              if(i>min(feat(find(a==k),2)))
                          t1= P(k,1)*i^2 + P(k,2)*i +P(k,3);
                          t2=round(t1);

                          %t3=max(t2-1,1);
                          %t4=min(t2+1,size(weight,1));
                          %weight(t3:t4,i)=weight(t3:t4,i)*penalty2; 


                                if(((t2-1)>=1)&&((t2+1)<size(weight,1)))
                                    weight(t2-1:t2+1,i)=penalty2;
                                end
                              end
                       end
                       end
                end
    %             matlabpool open 12;
                   parfor j=1:Q
                         Vn(:,j) = bp2(Phi, Xn(j,:).', weight(:,j),tau, opt.optNomalized);  
                         j
                   end

    %                matlabpool close;
    
    
