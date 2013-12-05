function mindex = multiguss(I,x,y)


 Ipart =I(x,y,:);
J = im2double(Ipart);


% for i = 1:a
%     for j = 1:b
%         sum = J(i,j,1)+J(i,j,2)+J(i,j,3);
%         if sum == 0
%             K(j,j,1) = 0;
%             K(j,j,2) = 0;
%             K(j,j,3) = 0;
%         else
%             K(i,j,1) = J(i,j,1)/sum;
%             K(i,j,2) = J(i,j,2)/sum;
%             K(i,j,3) = J(i,j,3)/sum;
%         end
%     end
% end


 R = J(:,:,1);
      G = J(:,:,2);
      B = J(:,:,3);

     [len,wid] = size(R);

     % Generation of 2-D Log Chromaticity Image.
     for i = 1:len
        for j = 1:wid
           if ((R(i,j)*G(i,j)*B(i,j))~= 0)
              c1(i,j) = R(i,j)/((R(i,j)*G(i,j)*B(i,j))^(1/3));
              c2(i,j) = G(i,j)/((R(i,j)*G(i,j)*B(i,j))^(1/3));
              c3(i,j) = B(i,j)/((R(i,j)*G(i,j)*B(i,j))^(1/3));
           else
              c1(i,j) = 1;
              c2(i,j) = 1;
              c3(i,j) = 1;
        end
    end
end

% rho1 = mat2gray(log(c1));
% rho2 = mat2gray(log(c2));
% rho3 = mat2gray(log(c3));
rho1 = (log(c1));
rho2 = (log(c2));
rho3 = (log(c3));

U1 = [1/sqrt(2),-1/sqrt(2),0];
U2 = [1/sqrt(6),1/sqrt(6),-2/sqrt(6)];

%  X1 = mat2gray((rho1-rho2)*(1/sqrt(2))); %(1/sqrt(2); -1/sqrt(2); 0)
%  X2 = mat2gray((rho1+rho2-2*rho3)*(1/(sqrt(6)))); %(1/sqrt(6); 1/sqrt(6); -2/sqrt(6))

 X1 = ((rho1-rho2)*(1/sqrt(2))); %(1/sqrt(2); -1/sqrt(2); 0)
 X2 = ((rho1+rho2-2*rho3)*(1/(sqrt(6)))); %(1/sqrt(6); 1/sqrt(6); -2/sqrt(6))

alltheta = 180;
%[h,p,jbstat,cv] = jbtest(x,alpha)siz
for t=1:alltheta
%     disp(t);
    delta = t*pi/alltheta;
    img = cos(delta)*X1 + sin(delta)*X2 ;
    %[h(t),p(t),jbstat(t),cv(t)] = lillietest(img(:,0.05,1));
    [h(t),p(t),jbstat(t),cv(t)] = jbtest(img(:),0.05,1);
    %hist(img);
%     logimg = -log(img);
%     r=img.*logimg; 
%     
%     n(t)=mean(abs(r(:)));
end;



%  plot(abs(h(:)), 'DisplayName', 'n', 'YDataSource', 'n'); figure(gcf)
%  figure ; plot(abs(p(:)), 'DisplayName', 'n', 'YDataSource', 'n');
%  figure; plot(X1,X2,'+');
 mindex = find(p == max(p));
 if length(mindex) >= 3
     plot(abs(p(:)), 'DisplayName', 'n', 'YDataSource', 'n');
 else
     mindex = mindex(1);
 end
 
 delta = mindex*pi/alltheta;
    img = cos(delta)*X1 + sin(delta)*X2 ;
%     figure;hist(img);

% for th = 50:10:130
%     L = cos(th*pi/180)*X1 + sin(th*pi/180)*X2;
%     imview(L);
% end


% theta = 145;
% L = cos(theta*pi/180)*X1 + sin(theta*pi/180)*X2;
% imtool(L);
% imwrite(L,'greyshadow.bmp');
% P = [cos(theta*pi/180),sin(theta*pi/180)]'*[cos(theta*pi/180),sin(theta*pi/180)];
% 
% [a,b] = size(X1);
% 
% 
% U = [U1',U2'];
% 
% for n1 = 1:a
%     for n2 = 1:b
%             C(n1,n2,:) = exp(U * P * [X1(n1,n2),X2(n1,n2)]');
%     end
% end
% 
% for n1 = 1:a
%     for n2 = 1:b
%         X(n1,n2,1) = C(n1,n2,1)/( C(n1,n2,1)+ C(n1,n2,2)+ C(n1,n2,3));
%         X(n1,n2,2) = C(n1,n2,2)/( C(n1,n2,1)+ C(n1,n2,2)+ C(n1,n2,3));
%         X(n1,n2,3) = C(n1,n2,3)/( C(n1,n2,1)+ C(n1,n2,2)+ C(n1,n2,3));
%     end
% end
% 
% %imview(X);
% imwrite(X,'colorshadow.jpg');
%       
% for n1 = 1:a
%     for n2 = 1:b
%         if (J(n1,n2,1)+ J(n1,n2,2)+ J(n1,n2,3)) == 0
%             LL(n1,n2,:) = 0;
%         else
%         LL(n1,n2,1) = J(n1,n2,1)/( J(n1,n2,1)+ J(n1,n2,2)+ J(n1,n2,3));
%         LL(n1,n2,2) = J(n1,n2,2)/( J(n1,n2,1)+ J(n1,n2,2)+ J(n1,n2,3));
%         LL(n1,n2,3) = J(n1,n2,3)/( J(n1,n2,1)+ J(n1,n2,2)+ J(n1,n2,3));
%         end
%     end
% end

        
        

% [v,idx]=min(n(:));