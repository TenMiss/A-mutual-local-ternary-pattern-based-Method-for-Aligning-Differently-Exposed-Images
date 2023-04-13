function [delta, phi] = TH_imfLTPalignment(im,TH1,TH2)
% Copyright: Wu Shiqian
% 11 Sep 2010
I = im{1};
[ri,ci,h] = size(I);
if h>1
    im{1} = rgb2gray(im{1});
    im{2} = rgb2gray(im{2});
end
A = [1,2,4,8,16,32,64,128,256,512];
B = ceil(min([ri,ci])/16);
D = B - A;
num = find(D>=0);
scale = 1./A(num);
n = length(num);
lp = fspecial('ga',5,1);
im0 = cell(1,n);
im1 = cell(1,n);
A = round(conv2(double(im{1}),lp,'same')); 
B = round(conv2(double(im{2}),lp,'same')); 
im0{1} = uint8(A);
im1{1} = uint8(B);
for i=2:n
    im0{i} = imresize(im0{1},scale(i),'bicubic');
    im1{i} = imresize(im1{1},scale(i),'bicubic');
end
%%%%%%%%%%%%%%%%%%%%%
gama =[];
for k = n:-1:n-2
    f0 = im0{k};
    f1 = im1{k};   
    [rr,cc]= size(f1);
    [ROTS,dxy,sigbeta] = TH_imfLTPCoarseTune(f0,f1,[-10,10],fix(rr/4),fix(cc/4),3,1,1,TH1,TH2);    
    if sigbeta==1   
        disp('Find the initial angle')
        break        
    else
        gama = [gama;ROTS];
    end
end

if sigbeta == 1
    nn = k-1;
    stot = [dxy(1),dxy(2),-mean(ROTS)];
    A1 = shift(im1{k-1},2*stot(2),2*stot(1)); 
    A1 = imrotate(A1,-stot(3),'bicubic','crop');
    srow = ceil(2*stot(1));
    scol = ceil(2*stot(2));
    A0 = im0{k-1};
    if srow > 0
        A0(1:srow,:)=[]; 
        A1(1:srow,:)=[]; 
    else
        A0(end+srow+1:end,:) =[];
        A1(end+srow+1:end,:) =[];        
    end
    if scol > 0
        A0(:,1:scol)=[];
        A1(:,1:scol)=[];
    else
        A0(:,end+scol+1:end) =[];
        A1(:,end+scol+1:end) =[];        
    end
    im0{k-1} = A0;
    im1{k-1} = A1;   
else
    nn = n-1;
    ROTS = median(gama(:));
    stot = [0,0,-mean(ROTS)];
    im1{n-1} = imrotate(im1{n-1},-stot(3),'bicubic','crop');
end   

for pyrlevel = nn:-1:1
    f0 = im0{pyrlevel};
    f1 = im1{pyrlevel}; 
    [f0, f1] = TH_imfLTPFeature(f0,f1,TH1,TH2);
    [y0,x0]=size(f0);
    xmean=x0/2; ymean = y0/2;
    x=kron((-xmean:xmean-1),ones(y0,1));
    y=kron(ones(1,x0),(-ymean:ymean-1)');
    sigma=1;
    g3 = exp(-(y.^2+x.^2)/(2*sigma^2))/2/pi/sigma^2;
    g1 = -g3.*y; 
    g2 = -g3.*x; 
    a=real(ifft2(fft2(f1).*fft2(g2))); 
    c=real(ifft2(fft2(f1).*fft2(g1))); 
    b=real(ifft2(fft2(f1).*fft2(g3)))-real(ifft2(fft2(f0).*fft2(g3))); 
    R=c.*x-a.*y; 
    a11 = sum(sum(a.*a)); a12 = sum(sum(a.*c)); a13 = sum(sum(R.*a));
    a21 = sum(sum(a.*c)); a22 = sum(sum(c.*c)); a23 = sum(sum(R.*c)); 
    a31 = sum(sum(R.*a)); a32 = sum(sum(R.*c)); a33 = sum(sum(R.*R));
    b1 = sum(sum(a.*b)); b2 = sum(sum(c.*b)); b3 = sum(sum(R.*b));
    
    Ainv = [a11 a12 a13; a21 a22 a23; a31 a32 a33];
    rk = rank(Ainv);
    if rk == 3
        Ainv = [a11 a12 a13; a21 a22 a23; a31 a32 a33]^(-1);
        s = Ainv*[b1; b2; b3];
        st = s;
        SGmax = max([abs(s(1)),abs(s(2)),abs(s(3))*180/pi]); 
        if SGmax <=100
            it=1;       
            while ((abs(s(1))+abs(s(2))+abs(s(3))*180/pi/20>0.01)&& it<=20)
                %while ((abs(s(1))+abs(s(2))+abs(s(3))*180/pi/20>0.002)&& it<20)
                if abs(st(3)*180/pi)>=1
                    tmp = shift(im0{pyrlevel},-st(1),-st(2));
                    tmp = imrotate(tmp,-st(3)*180/pi,'bicubic','crop');        
                    ff0 = TH_imfLTPFeature(tmp,im1{pyrlevel},TH1,TH2);     
                else
                    ff0 = shift(f0,-st(1),-st(2));
                end
                b = real(ifft2(fft2(f1).*fft2(g3)))-real(ifft2(fft2(ff0).*fft2(g3)));
                s = Ainv*[sum(sum(a.*b)); sum(sum(c.*b)); sum(sum(R.*b))]; 
                st = st+s;
                it = it+1;
                %err =abs(s(1))+abs(s(2))+abs(s(3))*180/pi/20;
            end
            st(3)=-st(3)*180/pi;
            st = st';
            st(1:2) = st(2:-1:1);
            stot = [2*stot(1:2)+st(1:2) stot(3)+st(3)];
            if pyrlevel>1
                im1{pyrlevel-1} = shift(im1{pyrlevel-1},2*stot(2),2*stot(1)); 
                im1{pyrlevel-1} = imrotate(im1{pyrlevel-1},-stot(3),'bicubic','crop');
                srow = ceil(2*stot(1));
                scol = ceil(2*stot(2));
                A0 = im0{pyrlevel-1};
                A1 = im1{pyrlevel-1};    
                if srow > 0
                    A0(1:srow,:)=[]; 
                    A1(1:srow,:)=[]; 
                else
                    A0(end+srow+1:end,:) =[];
                    A1(end+srow+1:end,:) =[];        
                end
                if scol > 0
                    A0(:,1:scol)=[];
                    A1(:,1:scol)=[];
                else
                    A0(:,end+scol+1:end) =[];
                    A1(:,end+scol+1:end) =[];        
                end
                im0{pyrlevel-1} = A0;
                im1{pyrlevel-1} = A1;   
            end
        end
    end
end
phi = stot(3);
delta = stot(1:2);
