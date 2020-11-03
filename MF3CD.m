clc; clear all; close all;

[fileName,FolderName] = uigetfile('*.*', 'Path selection Time 1');
cd(FolderName);
C = strsplit(FolderName,'_');
%date = cell2mat(C(8));
% date = 'r';

str = computer;

if str == 'PCWIN64'
    deli = '\';
elseif str == 'GLNXA64'
    deli = '/';
end
    

config_ID = fopen(strcat(FolderName,deli,'config.txt'),'rb');
tline = fgetl(config_ID);
tline = fgetl(config_ID);
b = str2num(tline); %row
% b = 936; %after reading from ENVI
tline = fgetl(config_ID);
tline = fgetl(config_ID);
tline = fgetl(config_ID);
a = str2num(tline); %column
% a = 979;%after reading from ENVI


cd(FolderName)
fileList = dir('*.bin');
%%
folderName = strcat(FolderName,'\','T11.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T11 = fread(fileID,[a b],'float32');
T11 = T11';

folderName = strcat(FolderName,'\','T12_imag.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T12_imag = fread(fileID,[a b],'float32');
T12_imag = T12_imag';

folderName = strcat(FolderName,'\','T12_real.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T12_real = fread(fileID,[a b],'float32');
T12_real = T12_real';

T12 = complex(T12_real,T12_imag);

folderName = strcat(FolderName,'\','T22.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T22 = fread(fileID,[a b],'float32');
T22 = T22';
%%
T21 = conj(T12);
%%
[nrows,ncols]= size(T11);
%%
% no optimization
disp('MF3CD')
db_theta = zeros(nrows,ncols);
pd = zeros(nrows,ncols);
pv = zeros(nrows,ncols);
ps = zeros(nrows,ncols);

% t11_avg = zeros(nrows,ncols);
% t22_avg = zeros(nrows,ncols);
% t33_avg = zeros(nrows,ncols);

dop_avg = zeros(nrows,ncols);
% wsi=input('Window Size: ');
wsi = 7;
disp('Taking window size as 7')
wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopj= nrows-inci; % Stop row for window processing
stopi= ncols-incj; % Stop column for window processing

t = cputime;

for ii=startj:stopj
    for jj=starti:stopi
        
        t11s = mean2(T11(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t12s = mean2(T12(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
                
        
        t21s = conj(mean2(T12(ii-inci:ii+inci,jj-incj:jj+incj)));%i sample
        t22s = mean2(T22(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        %Coherency matrix
        T_T1 = [t11s t12s; t21s t22s];
        
        m1 = real(sqrt(1-(4*(det(T_T1)./(trace(T_T1).^2))))); % DOP Barakat
        
        span = t11s + t22s;
        %         g1 = abs(t12s)^2+abs(t13s)^2;
        h = (t11s - t22s);
        %         h1 = (t11s + g1 - t22s - t33s);
        g = t22s;
        val = (m1.*span.*h)./(t11s.*g+m1.^2.*span.^2);
        
        eta = atand(val);
        db_theta(ii,jj) = eta;
        
        ps(ii,jj) = (((m1.*(span).*(1+sind(2*eta))./2)));
        pd(ii,jj) = (((m1.*(span).*(1-sind(2*eta))./2)));
        pv(ii,jj) = (span.*(1-m1));
        
        dop_avg(ii,jj) = m1;
    end
    fprintf('Column: %d \n',ii);
end
%%
path = FolderName;
if ~exist(strcat(path,'MF3CD'), 'dir')
    mkdir(strcat(path,'MF3CD'));
end
fclose('all');
FolderName = path;
Fold = strcat(FolderName,'MF3CD',deli);

path = Fold;
cd(path);
%%
addpath('E:\Prof_Juan_Opti_SAR_fusion\');
f_name_surface = strcat(['Ps_MF3CD','.bin']);
f_name_double_bounce = strcat(['Pd_MF3CD','.bin']);
f_name_diffused = strcat(['Pv_MF3CD','.bin']);

f_name_theta = strcat(['theta_MF3CD','.bin']);
f_name_dop = strcat(['dop_MF3CD','.bin']);

fileandpath_surface=strcat([path f_name_surface]);
fileandpath_double_bounce=strcat([path  f_name_double_bounce]);
fileandpath_diffused=strcat([path  f_name_diffused]);

fileandpath_theta=strcat([path  f_name_theta]);
fileandpath_dop=strcat([path  f_name_dop]);

fid_05 = fopen(fileandpath_surface,'wb');
fid_06 = fopen(fileandpath_double_bounce,'wb');
fid_07 = fopen(fileandpath_diffused,'wb');

fid_56 = fopen(fileandpath_theta,'wb');
fid_58 = fopen(fileandpath_dop,'wb');

fwrite(fid_05,ps', 'float32');
hdrwrite_envi('Ps_MF3CD', path, nrows, ncols)
fwrite(fid_06,pd', 'float32');
hdrwrite_envi('Pd_MF3CD', path, nrows, ncols)
fwrite(fid_07,pv', 'float32');
hdrwrite_envi('Pv_MF3CD', path, nrows, ncols)

fwrite(fid_56,db_theta', 'float32');
hdrwrite_envi('theta_MF3CD', path, nrows, ncols)
fwrite(fid_58,dop_avg', 'float32');
hdrwrite_envi('dop_MF3CD', path, nrows, ncols)

% close the file

fclose('all');

% imwrite(rgbVivid,'DB_Full_barakat.tif');
disp('End');