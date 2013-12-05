clear;
close all;
I = imread('1.jpg');
G = imcrop(I, [1,1,99,99]);%左上角与左下角
K = imcrop(I, [1,176,99,99]);

% I = imread('2.jpg');
% % G = imcrop(I, [1,1,99,99]);
% % K = imcrop(I, [253,189,99,99]);左上角与右下角
% 
% G = imcrop(I, [1,1,99,99]);
% K = imcrop(I, [1,189,99,99]);%左上角与左下角

Gd = im2double(G);
Kd = im2double(K);

[x,y,z] = size(Gd);

%化成log比值
for n1 = 1:x
    for n2 = 1:y
        div1 = Gd(n1,n2,1)*Gd(n1,n2,2)*Gd(n1,n2,3);
        div2 = Kd(n1,n2,1)*Kd(n1,n2,2)*Kd(n1,n2,3);
           if div1 ~= 0
        a(n1,n2) = log(Gd(n1,n2,2)/div1);
        b(n1,n2) = log(Gd(n1,n2,3)/div1);
           else
                a(n1,n2) = 1;
                b(n1,n2) = 1;
           end
           
            if div2 ~= 0
        c(n1,n2) = log(Kd(n1,n2,2)/div2);
        d(n1,n2) = log(Kd(n1,n2,3)/div2);
            else
        c(n1,n2) = 1;
        d(n1,n2) = 1;
           end
           
               
       % plot(a,b,'.');
   
    end
end


for th = 1:180
    posa = a*cos(th*pi/180) + b*sin(th*pi/180);%转成投影值
    posb = c*cos(th*pi/180) + d*sin(th*pi/180);

    [Xa,Ya] = hist(posa(:,:),-3:0.1:3);%60个区间的直方图
    [Xb,Yb] = hist(posb(:,:),-3:0.1:3);

    [s1,s2] = size(Xa);

    for i = 1:s1
        value(i) =  abs(sum(Xa(i,:)') - sum(Xb(i,:)'));%每个区间的直方图频数差
    end

    total(th) = sum(value);%每个角度的直方图相似度
end;

plot(abs(total(:)), 'DisplayName', 'n', 'YDataSource', 'n'); figure(gcf)
        
        