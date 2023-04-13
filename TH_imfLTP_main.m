% Image alignment + linear optimization
%
% Wu Shiqian. 11 Sep 2010

%% specify the directory that contains your range of differently exposed
clear all;
close all
clc
[filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Open images') %#ok<NOPTS>
if isequal([filename,pathname],[0,0]);
    return
end
thumbs_file = fullfile(pathname,'thumbs.db');
delete(thumbs_file)
[filenames, exposures, numExposures] = readDir(pathname);
% Pre-define parameters
tmp = imread(filenames{1});
I = cell(1,numExposures);
% Check the consistent of image labels
for i=1:numExposures
    a = imread(filenames{i});
    tmp = rgb2gray(a);
    I{i} = tmp;    
end

%% Image alignment
PM{1} = I{1};
[row,col,h]=size(PM);
shifts = zeros(numExposures-1,2);
alfa = zeros(1,numExposures-1);
TH1 = 2;
TH2 = 1;
PM{1} = I{1};
Ref = imread(filenames{1});
for i = 2:numExposures
    PM{2} = I{i};
    [sft, beta] = TH_imfLTPalignment(PM,TH1,TH2);
    shifts(i,:)=sft
    alfa(i) = beta;
    tmp = imread(filenames{i});
    if any(abs(shifts(i,:))>0.25)
        xx = floor(shifts(i,2)*4+0.5)/4;
        yy = floor(shifts(i,1)*4+0.5)/4;
        tmp = shift(tmp,xx,yy);
    end
    if abs(alfa(i))>0.01
        tmp = imrotate(tmp,-alfa(i),'bicubic','crop');
    end
%     name =['C:\Test\FamilyMod\' 'New' int2str(i) '.tiff'];
%     imwrite(tmp,name);
    figure, imshow(Ref); 
    hold on;
    h2 = imshow(tmp);
    title('overlay two aligned images'); 
    alpha(h2,0.6);
    hold off
end

alfa
shifts
fprintf('Finish image alignment by LTP\n')

Dalfa = abs(alfa + 5);
a = ones(numExposures,1)*[30,10];
Dshifts = abs(shifts-a);




