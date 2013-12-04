I = imread('1.jpg');
[a,b,c] = size(I);
step = 20;
adex = floor(a/step);
bdex = floor(b/step);

for i = 1:adex+1
    for j = 1:bdex
        if i == adex+1
            theta = multiguss(I,(i-1)*step+1:end,(j-1)*step+1:j*step);
        else
            
            theta = multiguss(I,(i-1)*step+1:i*step,(j-1)*step+1:j*step);
        end
        disp(theta);
    end

     if i == adex+1
        theta = multiguss(I,(i-1)*step+1:end,bdex*step+1:end);
     else
        theta = multiguss(I,(i-1)*step+1:i*step,bdex*step+1:end);
     end
    disp(theta);
end
        
