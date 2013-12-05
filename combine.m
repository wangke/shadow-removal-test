clc; clear all; close all;

rgb = imread('1.jpg');

H = fspecial('gaussian');
rgb = imfilter(rgb,H);

if ndims(rgb) == 3

    I = rgb2gray(rgb);

else

    I = rgb;

end

 

hy = fspecial('sobel');

hx = hy';

Iy = imfilter(double(I), hy, 'replicate');

Ix = imfilter(double(I), hx, 'replicate');

gradmag = sqrt(Ix.^2 + Iy.^2);

 

% L = watershed(gradmag);

% Lrgb = label2rgb(L);

 

se = strel('disk', 3);%此处是一个参数调整处

Io = imopen(I, se);

 

Ie = imerode(I, se);

Iobr = imreconstruct(Ie, I);

 

Ioc = imclose(Io, se);

 

Iobrd = imdilate(Iobr, se);

Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));

Iobrcbr = imcomplement(Iobrcbr);

 

fgm = imregionalmax(Iobrcbr);

 

It1 = rgb(:, :, 1);

It2 = rgb(:, :, 2);

It3 = rgb(:, :, 3);

It1(fgm) = 255; It2(fgm) = 0; It3(fgm) = 0;

I2 = cat(3, It1, It2, It3);

 

se2 = strel(ones(5,5));

fgm2 = imclose(fgm, se2);

fgm3 = imerode(fgm2, se2);

 

fgm4 = bwareaopen(fgm3, 20);

It1 = rgb(:, :, 1);

It2 = rgb(:, :, 2);

It3 = rgb(:, :, 3);

It1(fgm4) = 255; It2(fgm4) = 0; It3(fgm4) = 0;

I3 = cat(3, It1, It2, It3);

 

bw = im2bw(Iobrcbr, graythresh(Iobrcbr));

 

D = bwdist(bw);

DL = watershed(D);

bgm = DL == 0;

 

gradmag2 = imimposemin(gradmag, bgm | fgm4);

L = watershed(gradmag2);

 

It1 = rgb(:, :, 1);

It2 = rgb(:, :, 2);

It3 = rgb(:, :, 3);

fgm5 = imdilate(L == 0, ones(3, 3)) | bgm | fgm4;

It1(fgm5) = 255; It2(fgm5) = 0; It3(fgm5) = 0;

I4 = cat(3, It1, It2, It3);

 

Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
imtool(Lrgb);



for o = 1:max(L(:))  
    [x,y] = find(L == o);
    J = zeros(1,length(x),3);
    for p = 1:length(x)
        J(1,p,:) = rgb(x(p),y(p),:);
    end
    
  
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


    rho1 = (log(c1));
    rho2 = (log(c2));
    rho3 = (log(c3));

    X1 = ((rho1-rho2)*(1/sqrt(2))); %(1/sqrt(2); -1/sqrt(2); 0)
    X2 = ((rho1+rho2-2*rho3)*(1/(sqrt(6)))); %(1/sqrt(6); 1/sqrt(6); -2/sqrt(6))

    alltheta = 180;
    %[h,p,jbstat,cv] = jbtest(x,alpha)siz
    for t=1:alltheta
        disp(t);
        delta = t*pi/alltheta;
        img = cos(delta)*X1 + sin(delta)*X2 ;
        [h(t),p(t),kstat(t),cv(t)] = lillietest(img(:));
    end;


%     figure ; plot(abs(h(:)), 'DisplayName', 'n', 'YDataSource', 'n');
%     figure ; plot(abs(p(:)), 'DisplayName', 'n', 'YDataSource', 'n');
%     figure ; plot(abs(kstat(:)), 'DisplayName', 'n', 'YDataSource', 'n');
% 
%     figure; plot(X1,X2,'+');
    mindex(o) = find(kstat == min(kstat));
%     delta = mindex*pi/alltheta;
%     img = cos(delta)*X1 + sin(delta)*X2 ;
%     figure;hist(img);
end