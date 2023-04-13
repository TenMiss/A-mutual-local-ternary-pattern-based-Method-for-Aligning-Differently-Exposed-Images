function [IMFTable,Ref_His,Img_His] = IMF(Ref, Img)
% IMF mapping
% 5 Oct 2011
[row,col,h]=size(Ref);
if h > 1
    Ref = rgb2gray(Ref);
    Img = rgb2gray(Img);
end
Ref_His = zeros(1,256);
Img_His = zeros(1,256);
Ref = double(Ref(:))+1;
Img = double(Img(:))+1;
for i = 1:length(Ref)
    a = Ref(i);
    Ref_His(a) = Ref_His(a)+1;
    a = Img(i);
    Img_His(a) = Img_His(a)+1;
end
for i=2:256
    Ref_His(i) = Ref_His(i-1)+Ref_His(i);
    Img_His(i) = Img_His(i-1)+Img_His(i);
end
IMFTable = zeros(1,256);
for i=1:256
    a = abs(Ref_His(i)-Img_His);
    [v,num] = min(a);
    IMFTable(i)= num;    
end

