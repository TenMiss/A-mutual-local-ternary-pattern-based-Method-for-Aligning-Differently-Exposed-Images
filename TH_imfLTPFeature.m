function [FB, FD] = TH_imfLTPFeature(IB,ID,TH1,TH2)
% Extract local fetures FB, FD of images IB and ID
% Input: 
%        IB is a Bright Image with long exposure
%        ID is a Dark Image with short exposure
%        TH1 is a threshold for intensities(scalar)
%        TH2 is a threshold for IMF(scalar)
% Output:
%        FB,FD are local features of IB, ID with the same size of X
% Copyright: Wu Shiqian
% Institute for Infocomm Research
% 7 Oct 2011
% if nargin<3
%     TH1 = 2;     
%     TH2 = 2;   
% end
[row,col,h]=size(IB);
if h > 1
    IB = rgb2gray(IB);
    ID = rgb2gray(ID);
end
IB = double(IB); ID = double(ID);
f12 = IMF(IB,ID);
f21 = IMF(ID,IB);  
B = zeros(row+2,col+2); 
B2 = B; B4 = B; B6 = B; B8 = B;
B(2:row+1,2:col+1) = IB;
B2(3:row+2,2:col+1) = IB; 
B4(2:row+1,1:col) = IB; 
B6(1:row,2:col+1) = IB; 
B8(2:row+1,3:col+2) = IB; 
newB = f12(B+1)-1;
newB2 = f12(B2+1)-1;
newB4 = f12(B4+1)-1;
newB6 = f12(B6+1)-1;
newB8 = f12(B8+1)-1;

D = zeros(row+2,col+2); 
D2 = D; D4 = D; D6 = D; D8 = D;
D(2:row+1,2:col+1) = ID;
D2(3:row+2,2:col+1) = ID; 
D4(2:row+1,1:col) = ID; 
D6(1:row,2:col+1) = ID; 
D8(2:row+1,3:col+2) = ID; 
newD = f21(D+1)-1;
newD2 = f21(D2+1)-1;
newD4 = f21(D4+1)-1;
newD6 = f21(D6+1)-1;
newD8 = f21(D8+1)-1;

DA2 = B2 -B; DA4 = B4 -B; DA6 = B6 -B; DA8 = B8 -B; 
DB2 = newB2 -newB; DB4 = newB4 -newB; DB6 = newB6 -newB; DB8 = newB8 -newB; 
B1 = (DA2>TH1 & DB2>TH2); B2 = (DA2<-TH1 & DB2<-TH2);
BA2 = 2*B1+B2; 
B1 = (DA4>TH1 & DB4>TH2); B2 = (DA4<-TH1 & DB4<-TH2);
BA4 = 2*B1+B2; 
B1 = (DA6>TH1 & DB6>TH2); B2 = (DA6<-TH1 & DB6<-TH2);
BA6 = 2*B1+B2;
B1 = (DA8>TH1 & DB8>TH2); B2 = (DA8<-TH1 & DB8<-TH2);
BA8 = 2*B1+B2;
FB = BA2 + 3*BA4 + 9*BA6 + 27*BA8;
FB = FB(2:row+1,2:col+1);

DA2 = D2 -D; DA4 = D4 -D; DA6 = D6 -D; DA8 = D8 -D; 
DB2 = newD2 -newD; DB4 = newD4 -newD; DB6 = newD6 -newD; DB8 = newD8 -newD; 
B1 = (DA2>TH1 & DB2>TH2); B2 = (DA2<-TH1 & DB2<-TH2);
BA2 = 2*B1+B2; 
B1 = (DA4>TH1 & DB4>TH2); B2 = (DA4<-TH1 & DB4<-TH2);
BA4 = 2*B1+B2;
B1 = (DA6>TH1 & DB6>TH2); B2 = (DA6<-TH1 & DB6<-TH2);
BA6 = 2*B1+B2;
B1 = (DA8>TH1 & DB8>TH2); B2 = (DA8<-TH1 & DB8<-TH2);
BA8 = 2*B1+B2;
FD = BA2 + 3*BA4 + 9*BA6 + 27*BA8;
FD = FD(2:row+1,2:col+1);

% DA2 = B2 -B; DA4 = B4 -B; DA6 = B6 -B; DA8 = B8 -B; 
% DB2 = newB2 -newB; DB4 = newB4 -newB; DB6 = newB6 -newB; DB8 = newB8 -newB; 
% B1 = (DA2>TH1 & DB2>TH2); B2 = (abs(DA2)<TH1 & abs(DB2)<TH2);
% BA2 = 2*B1+B2; 
% B1 = (DA4>TH1 & DB4>TH2); B2 = (abs(DA4)<TH1 & abs(DB4)<TH2);
% BA4 = 2*B1+B2; 
% B1 = (DA6>TH1 & DB6>TH2); B2 = (abs(DA6)<TH1 & abs(DB6)<TH2);
% BA6 = 2*B1+B2;
% B1 = (DA8>TH1 & DB8>TH2); B2 = (abs(DA8)<TH1 & abs(DB8)<TH2);
% BA8 = 2*B1+B2;
% FB = BA2 + 3*BA4 + 9*BA6 + 27*BA8;
% FB = FB(2:row+1,2:col+1);
% 
% DA2 = D2 -D; DA4 = D4 -D; DA6 = D6 -D; DA8 = D8 -D; 
% DB2 = newD2 -newD; DB4 = newD4 -newD; DB6 = newD6 -newD; DB8 = newD8 -newD; 
% B1 = (DA2>TH1 & DB2>TH2); B2 = (abs(DA2)<TH1 & abs(DB2)<TH2);
% BA2 = 2*B1+B2; 
% B1 = (DA4>TH1 & DB4>TH2); B2 = (abs(DA4)<TH1 & abs(DB4)<TH2);
% BA4 = 2*B1+B2;
% B1 = (DA6>TH1 & DB6>TH2); B2 = (abs(DA6)<TH1 & abs(DB6)<TH2);
% BA6 = 2*B1+B2;
% B1 = (DA8>TH1 & DB8>TH2); B2 = (abs(DA8)<TH1 & abs(DB8)<TH2);
% BA8 = 2*B1+B2;
% FD = BA2 + 3*BA4 + 9*BA6 + 27*BA8;
% FD = FD(2:row+1,2:col+1);

