clc
clear all
close all
I = imread('1.jpg');  % ����ͼ��
I2 = rgb2gray(I);
imshow(I);title('ԭͼ')
BW1 = edge(I2,'canny');  % ����canny����
figure,imshow(BW1);  % ��ʾ�ָ���ͼ�񣬼��ݶ�ͼ��
title('Canny')

 se = strel('square', 5);
 BW2 = imclose(BW1, se);%�ղ���
 figure,imshow(BW2);
% 
% BW2 = 1-BW2;
% %BW2 = im2bw(BW2);
% [BW3,num] = bwlabel(BW2);%��ͨ�Է���
% 
% 
% 
% %imtool(BW3>1);
% hold on;
% y = size(I,1);
% x = size(I,2);
% result = 0;
% index = 0;
% stats = regionprops(BW3, 'Centroid');
% for i = 1 : length(stats)
%     temp = stats(i).Centroid;
%     plot(temp(1), temp(2), 'r.');
%     way = ((temp(1)-x/2).^2 + (temp(2)-y/2).^2).^(1/2);
%     changeWay = 1/(1+exp(way*(-1)));%�����˹��Ȩ
%     area = sum(sum(BW3 == i));
%     if i == 1
%         result = changeWay*area;
%         index = i;
%     else
%         if changeWay*area > result
%             result = changeWay*area;
%             index = i;
%         end
%     end
% end
% 
% imtool(BW3 == index);
% 
% [x_s, y_s] = find(BW3 == index);
% 
% R_mean = mean(diag(I(x_s(:),y_s(:),1)));
% G_mean = mean(diag(I(x_s(:),y_s(:),2)));
% B_mean = mean(diag(I(x_s(:),y_s(:),3)));
% 
% R_var = var(double(diag(I(x_s(:),y_s(:),1))));
% G_var = var(double(diag(I(x_s(:),y_s(:),2))));
% B_var = var(double(diag(I(x_s(:),y_s(:),3))));
% 
% 
% 
% k = -1;
% 
% colordoor = 30;
% vardoor = 20;
% 
% for i = 1:num
%     
%         disp(i);
%         [x_s, y_s] = find(BW3 == i);
%            r_mean = mean(diag(I(x_s(:),y_s(:),1)));
%            g_mean = mean(diag(I(x_s(:),y_s(:),2)));
%             b_mean = mean(diag(I(x_s(:),y_s(:),3)));
% 
% r_var = var(double(diag(I(x_s(:),y_s(:),1))));
% g_var = var(double(diag(I(x_s(:),y_s(:),2))));
% b_var = var(double(diag(I(x_s(:),y_s(:),3))));
% 
%            %disp(r_var);disp(g_var);disp(b_var);
%            if abs(r_mean - R_mean) < colordoor && abs(g_mean - G_mean) < colordoor && abs(r_mean - R_mean) < colordoor
%                disp('ok');
%                if k == -1
%                    k = i;
%                else
%                    k = [k,i];
%                end
%            end
% 
% end
% 
% a = zeros(y,x);
% for i = 1:length(k)
%     a = a + (BW3 == k(i));
% end
% 
% %a = a >=1;
% imtool(a);


