clc;
clear all;
close all;

a = 0;
b = 0;

[filename, path] = uigetfile('*.txt', 'Select Full Pol config file');
%     FolderName = 'D:\FP_DP_CP_comp_RS2\processed\RS2_OK106814_PK930269_DK875264_FQ15W_20190606_124518_HH_VV_HV_VH_SLC_LEE\T3';
FolderName = path;
cd(FolderName);
% path = FolderName;
C = strsplit(FolderName,'\');
words = numel(C);
string = C(1,1);
save_dir = strcat(string,'\');
for i = 2:(words-2)
    str = C(1,i);
    save_dir = strcat(save_dir,str,'\');
end
fprintf('Save directory: %s\n', char(save_dir))
%     str = C(1);
%     save_dir = strcat(str,'\');
%     save_dir = string(save_dir);

% date = cell2mat(C(4));
date = 'r';

config_ID = fopen(strcat(FolderName,'\','config.txt'),'rb');
tline = fgetl(config_ID);
tline = fgetl(config_ID);
b = str2num(tline); %row
tline = fgetl(config_ID);
tline = fgetl(config_ID);
tline = fgetl(config_ID);
a = str2num(tline); %column
nrow = b;
ncol = a;

%     cd(FolderName)
fileList = dir('*.bin');

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
T21 = conj(T12);

folderName = strcat(FolderName,'\','T13_imag.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T13_imag = fread(fileID,[a b],'float32');
T13_imag = T13_imag';

folderName = strcat(FolderName,'\','T13_real.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T13_real = fread(fileID,[a b],'float32');
T13_real = T13_real';

T13 = complex(T13_real,T13_imag);
T31 = conj(T13);

folderName = strcat(FolderName,'\','T22.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T22 = fread(fileID,[a b],'float32');
T22 = T22';

folderName = strcat(FolderName,'\','T23_imag.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T23_imag = fread(fileID,[a b],'float32');
T23_imag = T23_imag';

folderName = strcat(FolderName,'\','T23_real.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T23_real = fread(fileID,[a b],'float32');
T23_real = T23_real';

T23 = complex(T23_real,T23_imag);
T32 = conj(T23);

folderName = strcat(FolderName,'\','T33.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T33 = fread(fileID,[a b],'float32');
T33 = T33';

% pd_f = zeros(nrow, ncol); % full
% pv_f = zeros(nrow, ncol); % full
% ps_f = zeros(nrow, ncol); % full
% theta_val_f = zeros(nrow, ncol); % full
% dop_val_f = zeros(nrow, ncol); % full

disp('load complete');

% pd_d = zeros(nrow, ncol); % dual co
% pv_d = zeros(nrow, ncol); % dual co
% ps_d = zeros(nrow, ncol); % dual co
% theta_val_d = zeros(nrow, ncol); % dual co
% dop_val_d = zeros(nrow, ncol); % dual co

wsi = 7;
disp('Taking window size as 7')

kernel = ones(wsi, wsi)./(wsi*wsi);
t11s = conv2(T11,kernel,'same');
t12s = conv2(T12,kernel,'same');
t13s = conv2(T13,kernel,'same');
t21s = conv2(T21,kernel,'same');
t22s = conv2(T22,kernel,'same');
t23s = conv2(T23,kernel,'same');
t31s = conv2(T31,kernel,'same');
t32s = conv2(T32,kernel,'same');
t33s = conv2(T33,kernel,'same');

T2_det = (t11s.*t22s-t12s.*t21s);
T2_trace = t11s+t22s;
s0_d = T2_trace;

dop_val_d = real(sqrt(1.0-(4.0.*T2_det./(T2_trace.^2))));

h = (t11s - t22s);
g = t11s;
val_d = (dop_val_d.*T2_trace.*h)./((t22s).*g+dop_val_d.^2.*T2_trace.^2);

theta_val_d = atand(val_d);
pd_d = dop_d.*(s0_d./2).*(1-sind(2*theta_d));
pv_d = (1-dop_d).*s0_d;
ps_d = dop_d.*(s0_d./2).*(1+sind(2*theta_d));

% wsj = wsi; % Number of columns in the window
% 
% inci=fix(wsi/2); % Up & down movement margin from the central row
% incj=fix(wsj/2); % Left & right movement from the central column
% % Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999
% 
% starti=fix(wsi/2)+1; % Starting row for window processing
% startj=fix(wsj/2)+1; % Starting column for window processing
% 
% stopj= nrow-inci; % Stop row for window processing
% stopi= ncol-incj; % Stop column for window processing
% 
% for ii = startj:stopj
%   for jj = starti:stopi
%       t11s_d = T11(ii-inci:ii+inci,jj-incj:jj+incj);
%       t11s_d = mean(t11s_d(:));
%       t12s_d = T12(ii-inci:ii+inci,jj-incj:jj+incj);
%       t12s_d = mean(t12s_d(:));
%       t21s_d = T21(ii-inci:ii+inci,jj-incj:jj+incj);
%       t21s_d = mean(t21s_d(:));
%       t22s_d = T22(ii-inci:ii+inci,jj-incj:jj+incj);
%       t22s_d = mean(t22s_d(:));
%       
%       T_d = [t11s_d, t12s_d; t21s_d, t22s_d];
%       
%       s0_d = trace(T_d);
%       
%       dop_d = real(sqrt(1-(4*(det(T_d)./(trace(T_d).^2)))));
%       dop_val_d(ii,jj) = dop_d;
% 
%       span = t11s_d + t22s_d;
%       h = (t11s_d - t22s_d);
%       g = t11s_d;
%       val_d = (dop_d.*span.*h)./((t22s_d).*g+dop_d.^2.*span.^2);
%       
%       theta_d = atand(val_d);
%       theta_val_d(ii,jj) = theta_d;
%       
%       pd_d(ii,jj) = dop_d.*(s0_d./2).*(1-sind(2*theta_d));
%       pv_d(ii,jj) = (1-dop_d).*s0_d;
%       ps_d(ii,jj) = dop_d.*(s0_d./2).*(1+sind(2*theta_d));
%   end
%   fprintf('row: %d\n',ii);
% end
%%
% addpath('D:\11. FP-DP-CP_comp_rice_RSE\');
addpath('ADD hdrwrite_envi.m path here');
fName = strcat('DoP_DP', '_', char(date));
f_name_dopDP = strcat('\', char(fName),'.bin');
fileandpath_dopDP=strcat([path f_name_dopDP]);
fid_021 = fopen(fileandpath_dopDP,'wb');
fwrite(fid_021,dop_val_d', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('DoP_DP bin and hdr files written');

fName = strcat('theta_DP', '_', char(date));
f_name_thetaDP = strcat('\', char(fName), '.bin');
fileandpath_thetaDP=strcat([path f_name_thetaDP]);
fid_02 = fopen(fileandpath_thetaDP,'wb');
fwrite(fid_02,theta_val_d', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('theta_DP bin and hdr files written');

fName = strcat('pd_DP', '_', char(date));
f_name_pdDP = strcat('\', char(fName), '.bin');
fileandpath_pdDP=strcat([path f_name_pdDP]);
fid_05 = fopen(fileandpath_pdDP,'wb');
fwrite(fid_05,pd_d', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('pd_DP bin and hdr files written');

fName = strcat('pv_DP', '_', char(date));
f_name_pvDP = strcat('\', char(fName), '.bin');
fileandpath_pvDP=strcat([path f_name_pvDP]);
fid_08 = fopen(fileandpath_pvDP,'wb');
fwrite(fid_08,pv_d', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('pv_DP bin and hdr files written');

fName = strcat('ps_DP', '_', char(date));
f_name_psDP = strcat('\', char(fName), '.bin');
fileandpath_psDP=strcat([path f_name_psDP]);
fid_11 = fopen(fileandpath_psDP,'wb');
fwrite(fid_11,ps_d', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('ps_DP bin and hdr files written');

fclose('all');
%%
cd(path)
fl1 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
set(fl1,'name','Theta','numbertitle','off');
imagesc(theta_val_d);
axis('image');
axis off;
caxis([-45 45]);
% title('Theta');
colormap(flipud(jet));
colorbar('FontSize', 20);

file1 =  char(strcat(path,'\theta_DP', '_', char(date)));
saveas(fl1,file1,'png')
saveas(fl1,file1,'fig')

disp('theta_DP.png saved');
