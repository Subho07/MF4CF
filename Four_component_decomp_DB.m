[filename, path] = uigetfile('*.txt', 'Select Full Pol config file');

FolderName = path;
%     cd(FolderName);
% path = FolderName;

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

pd_f = zeros(nrow, ncol); % full
pv_f = zeros(nrow, ncol); % full
ps_f = zeros(nrow, ncol); % full
pc_f = zeros(nrow, ncol); % full
theta_val_f = zeros(nrow, ncol); % full
tau_val_f = zeros(nrow, ncol); % full
dop_val_f = zeros(nrow, ncol); % full

disp('load complete');

%%
wsi = 7;
disp('Taking window size as 7')
wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopj= nrow-inci; % Stop row for window processing
stopi= ncol-incj; % Stop column for window processing

%%
for ii = startj:stopj
  for jj = starti:stopi
    
    
        t11s = T11(ii-inci:ii+inci,jj-incj:jj+incj);
        t11s = mean(t11s(:));
        t12s = T12(ii-inci:ii+inci,jj-incj:jj+incj);
        t12s = mean(t12s(:));
        t21s = T21(ii-inci:ii+inci,jj-incj:jj+incj);
        t21s = mean(t21s(:));
        t13s = T13(ii-inci:ii+inci,jj-incj:jj+incj);
        t13s = mean(t13s(:));
        t31s = T31(ii-inci:ii+inci,jj-incj:jj+incj);
        t31s = mean(t31s(:));
        t22s = T22(ii-inci:ii+inci,jj-incj:jj+incj);
        t22s = mean(t22s(:));
        t23s = T23(ii-inci:ii+inci,jj-incj:jj+incj);
        t23s = mean(t23s(:));
        t32s = T32(ii-inci:ii+inci,jj-incj:jj+incj);
        t32s = mean(t32s(:));
        t33s = T33(ii-inci:ii+inci,jj-incj:jj+incj);
        t33s = mean(t33s(:));

        T = [t11s, t12s, t13s; t21s, t22s, t23s; t31s, t32s, t33s];

        dop_f = real(sqrt(1-(27*(det(T)./(trace(T).^3)))));
        dop_val_f(ii,jj) = dop_f;

        k11_f = (t11s + t22s + t33s)/2;
        k44_f = (-t11s + t22s + t33s)/2;
        k14_f = imag(t23s);
%         k14_f = (-1i*(t23s - t32s))./2;
        
%         span = t11s + t22s + t33s;
%         %         g1 = abs(t12s)^2+abs(t13s)^2;
%         h = (t22s + t33s - k14_f);
%         %         h1 = (t11s + g1 - t22s - t33s);
%         g = k14_f;
%         val = (dop_f.*span.*h)./((t22s + t33s).*g+dop_f.^2.*span.^2);
%         
% %         eta = atand(val);
        

        s0_f = trace(T);

        
        val1 = (4*dop_f*k11_f*k44_f)./(k44_f^2 - (1 + 4*dop_f^2)*k11_f^2);
        val2 = abs(k14_f)./(k11_f);
       
        theta_f = atand(val1); % separation for surface and dbl
        tau_f = atand(val2); % separation for helix
        theta_val_f(ii,jj) = theta_f;
        tau_val_f(ii,jj) = tau_f;
        
        pv_f(ii,jj) = (1-dop_f).*s0_f;
        ps_f(ii,jj) = dop_f.*(s0_f./2).*(1+sind(2*theta_f));
        
        if theta_f < 0
            pc_f(ii,jj) = dop_f.*s0_f.*sind(2.*tau_f);
            pd_f(ii,jj) = s0_f - (pv_f(ii,jj) + ps_f(ii,jj) + pc_f(ii,jj));
        elseif theta_f >= 0
            theta_f_dbl = theta_f + tau_f;
            pd_f(ii,jj) = dop_f.*(s0_f./2).*(1-sind(2*theta_f_dbl));
            pc_f(ii,jj) = s0_f - (pv_f(ii,jj) + ps_f(ii,jj) + pd_f(ii,jj));
        end       
  end
  fprintf('row: %d\n',ii);
end
%%
path = FolderName;
if ~exist(strcat(path,'Four_component_DB_3'), 'dir')
    mkdir(strcat(path,'Four_component_DB_3'));
end
fclose('all');
FolderName = path;
Fold = strcat(FolderName,'Four_component_DB_3\');

path = Fold;
cd(path);
%%
% dop
% fig1 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% gamma100 = imagesc(dop_val_f);
% set(gamma100,'AlphaData',~isnan(dop_val_f))
% axis('image');
% axis off;
% 
% c = jet;
% colormap(c);
% colorbar('YTick', 0:0.2:1,'FontSize', 30);
% caxis([0 1]);

fig1 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
set(fig1,'name','DoP','numbertitle','off');
imagesc(dop_val_f);
axis('image');
axis off;
caxis([0 1]);
% title('Lambda 1');
colormap(jet);
colorbar('FontSize', 20);

file1 =  char('dop_Full_barakat');
saveas(fig1,file1,'png')
saveas(fig1,file1,'fig')

% theta
% fig2 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% gamma100 = imagesc(theta_val_f);
% set(gamma100,'AlphaData',~isnan(theta_val_f))
% axis('image');
% axis off;
% 
% c = jet;
% c = flipud(c);
% colormap(c);
% colorbar('YTick', -45:15:45,'FontSize', 30);
% caxis([-45 45]);

fig2 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
set(fig2,'name','Theta','numbertitle','off');
imagesc(theta_val_f);
axis('image');
axis off;
caxis([-45 45]);
% title('Lambda 1');
colormap(flipud(jet));
colorbar('FontSize', 20)

file1 =  char('theta_Full_barakat');
saveas(fig2,file1,'png')
saveas(fig2,file1,'fig')

% tau
% fig3 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% tau_val_f1 = tau_val_f;
% % tau_val_f1 = imadjust(tau_val_f,stretchlim(tau_val_f),[]);
% % min_val = min(tau_val_f1(:));
% min_val = 0;
% % max_val = max(tau_val_f1(:));
% max_val = 45;
% gamma100 = imshow(tau_val_f1);
% set(gamma100,'AlphaData',~isnan(tau_val_f1))
% axis('image');
% axis off;

% c = jet;
% % c = flipud(c);
% colormap(c);
% colorbar('YTick', min_val:((max_val-min_val)/5):max_val,'FontSize', 30);
% caxis([min_val max_val]);

fig3 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
tau_val_f1 = tau_val_f;
set(fig3,'name','Tau','numbertitle','off');
imagesc(tau_val_f1);
axis('image');
axis off;
caxis([0 45]);
% title('Lambda 1');
colormap(jet);
colorbar('FontSize', 20)

file1 =  char('tau_Full_barakat');
saveas(fig3,file1,'png')
saveas(fig3,file1,'fig')
%%
addpath('D:\27. Four_component_decomp_DB\');
fName = 'DoP_FP';
f_name_dopFP = strcat('\', char(fName), '.bin');
fileandpath_dopFP=strcat([path f_name_dopFP]);
fid_011 = fopen(fileandpath_dopFP,'wb');
fwrite(fid_011,dop_val_f', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('DoP_FP bin and hdr files written')

fName = 'theta_FP';
f_name_thetaFP = strcat('\', char(fName), '.bin');
fileandpath_thetaFP=strcat([path f_name_thetaFP]);
fid_01 = fopen(fileandpath_thetaFP,'wb');
fwrite(fid_01,theta_val_f', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('theta_FP bin and hdr files written');

fName = 'tau_FP';
f_name_tauFP = strcat('\', char(fName), '.bin');
fileandpath_tauFP=strcat([path f_name_tauFP]);
fid_01 = fopen(fileandpath_tauFP,'wb');
fwrite(fid_01,tau_val_f', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('tauFP bin and hdr files written');

fName = 'pd_FP';
f_name_pdFP = strcat('\', char(fName), '.bin');
fileandpath_pdFP=strcat([path f_name_pdFP]);
fid_04 = fopen(fileandpath_pdFP,'wb');
fwrite(fid_04,pd_f', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('pd_FP bin and hdr files written');

fName = 'pv_FP';
f_name_pvFP = strcat('\', char(fName), '.bin');
fileandpath_pvFP=strcat([path f_name_pvFP]);
fid_07 = fopen(fileandpath_pvFP,'wb');
fwrite(fid_07,pv_f', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('pv_FP bin and hdr files written');

fName = 'ps_FP';
f_name_psFP = strcat('\', char(fName), '.bin');
fileandpath_psFP=strcat([path f_name_psFP]);
fid_10 = fopen(fileandpath_psFP,'wb');
fwrite(fid_10,ps_f', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('ps_FP bin and hdr files written');

fName = 'pc_FP';
f_name_pcFP = strcat('\', char(fName), '.bin');
fileandpath_pcFP=strcat([path f_name_pcFP]);
fid_11 = fopen(fileandpath_pcFP,'wb');
fwrite(fid_11,pc_f', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('pc_FP bin and hdr files written');

im1 = cat(4,pd_f',pv_f',ps_f',pc_f');
f_name_allpowers = strcat(['Pd_Pv_Ps_Pc','.bin']);
fileandpath_allpowers = strcat([path f_name_allpowers]);
fid_09 = fopen(fileandpath_allpowers,'wb');
fwrite(fid_09,im1, 'float32');

fclose('all');

%%
warning('off');
clc;
close all;
fontSize = 16;
imshow(theta_val_f, []);
axis on;
title('Original Grayscale Image', 'FontSize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

% message = sprintf('Left click to select top pixel');
% uiwait(msgbox(message));
hFH = imfreehand();
% Create a binary image ("mask") from the ROI object.
% binaryImage = hFH.createMask();
xy = hFH.getPosition;
top_pos = xy(1,:);
% message = sprintf('Left click to select bottom pixel');
% uiwait(msgbox(message));
hFH = imfreehand();
% Create a binary image ("mask") from the ROI object.
% binaryImage1 = hFH.createMask();
xy = hFH.getPosition;
bottom_pos = xy(1,:);
% ROI x and y co-ordinates
top_x = top_pos(1,1);
top_y = top_pos(1,2);
top_pos = int32([top_x top_y]);
top_y = top_pos(1,1);
top_x = top_pos(1,2);
% h = drawcircle('Center',[top_x,top_y],'Radius', 6,'StripeColor','red');

bottom_x = bottom_pos(1,1);
bottom_y = bottom_pos(1,2);
bottom_pos = int32([bottom_x bottom_y]);
bottom_y = bottom_pos(1,1);
bottom_x = bottom_pos(1,2);

% h1 = drawcircle('Center',[bottom_x,bottom_y],'Radius', 6,'StripeColor','red');
% subsampled as per ROI
pd_f_sub = pd_f(top_x:bottom_x,top_y:bottom_y);
ps_f_sub = ps_f(top_x:bottom_x,top_y:bottom_y);
pc_f_sub = pc_f(top_x:bottom_x,top_y:bottom_y);
pv_f_sub = pv_f(top_x:bottom_x,top_y:bottom_y);
dop_val_f_sub = dop_val_f(top_x:bottom_x,top_y:bottom_y);
theta_val_f_sub = theta_val_f(top_x:bottom_x,top_y:bottom_y);
tau_val_f_sub = tau_val_f(top_x:bottom_x,top_y:bottom_y);

% average
pd_f_mean = mean(pd_f_sub(:));
ps_f_mean = mean(ps_f_sub(:));
pc_f_mean = mean(pc_f_sub(:));
pv_f_mean = mean(pv_f_sub(:));
dop_val_f_mean = mean(dop_val_f_sub(:));
theta_val_f_mean = mean(theta_val_f_sub(:));
tau_val_f_mean = mean(tau_val_f_sub(:));

% percentages
pd_f_percent = pd_f_mean./(pd_f_mean + ps_f_mean + pc_f_mean + pv_f_mean);
ps_f_percent = ps_f_mean./(pd_f_mean + ps_f_mean + pc_f_mean + pv_f_mean);
pc_f_percent = pc_f_mean./(pd_f_mean + ps_f_mean + pc_f_mean + pv_f_mean);
pv_f_percent = pv_f_mean./(pd_f_mean + ps_f_mean + pc_f_mean + pv_f_mean);

X = [pd_f_percent,ps_f_percent,pc_f_percent,pv_f_percent];

fig3 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
labels = {' ',' ',' ',' '};
pie(X,labels);
colormap([1 0 0;      %// red
          0 0 1;      %// blue
          1 0.08 0.6;      %// red_deriv
          0 1 0])  %// green

user_entry = input(' give file name: ', 's');
file1 =  char(user_entry);
saveas(fig3,file1,'png')
saveas(fig3,file1,'fig')

save_txt = strcat(user_entry,'.txt');
fid = fopen(save_txt,'wt');

for jj = 1:11
    if jj == 1
        str_p = 'pd';
    end
    if jj == 2
        str_p = 'ps';
    end
    if jj == 3
        str_p = 'pc';
    end
    if jj == 4
        str_p = 'pv';
    end
    if jj == 5
        str_p = 'dop';
    end
    if jj == 6
        str_p = 'theta';
    end
    if jj == 7
        str_p = 'tau';
    end
    if jj == 8
        str_p = 'top_x';
    end
    if jj == 9
        str_p = 'top_y';
    end
    if jj == 10
        str_p = 'bottom_x';
    end
    if jj == 11
        str_p = 'bottom_y';
    end
    fprintf(fid,'%s\t',str_p);
end

fprintf(fid,'\n');

for ii = 1:size(X,1)
    fprintf(fid,'%g\t',X(ii,:).*100);
%     fprintf(fid,'\n');
end

for ii = 1:3
    if ii == 1
        fprintf(fid,'%g\t',dop_val_f_mean);
    end
    if ii == 2
        fprintf(fid,'%g\t',theta_val_f_mean);
    end
    if ii == 3
        fprintf(fid,'%g\t',tau_val_f_mean);
    end
end

for ii = 1:4
    if ii == 1
        fprintf(fid,'%g\t',top_x);
    end
    if ii == 2
        fprintf(fid,'%g\t',top_y);
    end
    if ii == 3
        fprintf(fid,'%g\t',bottom_x);
    end
    if ii == 4
        fprintf(fid,'%g\t',bottom_y);
    end
%     fprintf(fid,'\n');
end
fprintf(fid,'\n');
fclose(fid);

% dom_class_DB_1 = zeros(nrow,ncol);
% dom_class_DB_2 = zeros(nrow,ncol);
% dom_class_DB_3 = zeros(nrow,ncol);
% dom_class_DB_4 = zeros(nrow,ncol);
% % for ii = 1:nrow
% %     for jj = 1:ncol
% %         if max([pd_f(ii,jj),ps_f(ii,jj),pv_f(ii,jj),pc_f(ii,jj)])== pd_f(ii,jj)
% %             dom_class_DB(ii,jj) = 1;
% %         end
% %         if max([pd_f(ii,jj),ps_f(ii,jj),pv_f(ii,jj),pc_f(ii,jj)])== ps_f(ii,jj)
% %             dom_class_DB(ii,jj) = 2;
% %         end
% %         if max([pd_f(ii,jj),ps_f(ii,jj),pv_f(ii,jj),pc_f(ii,jj)])== pv_f(ii,jj)
% %             dom_class_DB(ii,jj) = 3;
% %         end
% %         if max([pd_f(ii,jj),ps_f(ii,jj),pv_f(ii,jj),pc_f(ii,jj)])== pc_f(ii,jj)
% %             dom_class_DB(ii,jj) = 4;
% %         end
% %     end
% % end
% 
% for ii = 1:nrow
%     for jj = 1:ncol
%         vec = [pd_f(ii,jj),ps_f(ii,jj),pv_f(ii,jj),pc_f(ii,jj)];
%         [B,I] = sort(vec,'descend');
% %         dom_class_DB_1(ii,jj) = I(1,1);
% %         dom_class_DB_2(ii,jj) = I(1,2);
% %         dom_class_DB_3(ii,jj) = I(1,3);
% %         dom_class_DB_4(ii,jj) = I(1,4);
%         dom = I(1,1);
%         dom_norm = B(1,1)/sum(B);
%         if dom == 1
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_DB_1(ii,jj) = 1;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_DB_1(ii,jj) = 2;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_DB_1(ii,jj) = 3;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_DB_1(ii,jj) = 4;
%             end
%         end
%         
%         if dom == 2
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_DB_1(ii,jj) = 5;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_DB_1(ii,jj) = 6;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_DB_1(ii,jj) = 7;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_DB_1(ii,jj) = 8;
%             end
%         end
%         
%         if dom == 3
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_DB_1(ii,jj) = 9;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_DB_1(ii,jj) = 10;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_DB_1(ii,jj) = 11;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_DB_1(ii,jj) = 12;
%             end
%         end
%         
%         if dom == 4
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_DB_1(ii,jj) = 13;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_DB_1(ii,jj) = 14;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_DB_1(ii,jj) = 15;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_DB_1(ii,jj) = 16;
%             end
%         end
%             
%     end
% end

% f7 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% set(gca,'FontSize',20)
% imagesc(dom_class_DB_1)
% axis('image');
% axis off;
% 
% mymap = [0.4,0,0;
%     1,0,0;
%     0.88,0.24,0.16;
%     0.80,0.36,0.36;%red
%     0.53,0.1,0.98;
%     0.12,0.56,1;
%     0.48,0.41,0.93;
%     0,0,1;%blue
%     0.68,1,0.18;
%     0.56,0.74,0.56;
%     0.13,0.54,0.13;
%     0,1,0;%green
%     1,0.08,0.57;
%     1,0.41,0.71;
%     1,0.71,0.76;
%     1,0.74,0.79;%pink
%    ];
% 
% colormap(mymap);
% 
% numcolors = 16;
% caxis([1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% cbarHandle = colorbar('YTick',...
% [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% 'YTickLabel',{'Z1','Z2','Z3','Z4','Z5','Z6','Z7','Z8','Z9','Z10','Z11','Z12','Z13','Z14','Z15','Z16'}, 'YLim', [1 numcolors]);
% set(gca,'FontSize', 20);
% file1 =  char(strcat(path,'\DB_dom_1'));
% saveas(f7,file1,'png')
% saveas(f7,file1,'fig')
% 
% % f8 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% % set(gca,'FontSize',20)
% % imagesc(dom_class_DB_2)
% % axis('image');
% % axis off;
% % 
% % mymap = [1,0,0;
% %     0,0,1;
% %     0,1,0;
% %     1,0.76,0.8;
% %    ];
% % 
% % colormap(mymap);
% % 
% % numcolors = 4;
% % caxis([1 numcolors]);
% % % cbarHandle = colorbar('YTick',...
% % % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',{'Even','Odd','Diff','Hlx'}, 'YLim', [1 numcolors]);
% % set(gca,'FontSize', 20);
% % file2 =  char(strcat(path,'\DB_dom_2'));
% % saveas(f8,file2,'png')
% % saveas(f8,file2,'fig')
% % 
% % f9 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% % set(gca,'FontSize',20)
% % imagesc(dom_class_DB_3)
% % axis('image');
% % axis off;
% % 
% % mymap = [1,0,0;
% %     0,0,1;
% %     0,1,0;
% %     1,0.76,0.8;
% %    ];
% % 
% % colormap(mymap);
% % 
% % numcolors = 4;
% % caxis([1 numcolors]);
% % % cbarHandle = colorbar('YTick',...
% % % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',{'Even','Odd','Diff','Hlx'}, 'YLim', [1 numcolors]);
% % set(gca,'FontSize', 20);
% % file3 =  char(strcat(path,'\DB_dom_3'));
% % saveas(f9,file3,'png')
% % saveas(f9,file3,'fig')
% % 
% % f10 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% % set(gca,'FontSize',20)
% % imagesc(dom_class_DB_4)
% % axis('image');
% % axis off;
% % 
% % mymap = [1,0,0;
% %     0,0,1;
% %     0,1,0;
% %     1,0.76,0.8;
% %    ];
% % 
% % colormap(mymap);
% % 
% % numcolors = 4;
% % caxis([1 numcolors]);
% % % cbarHandle = colorbar('YTick',...
% % % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',{'Even','Odd','Diff','Hlx'}, 'YLim', [1 numcolors]);
% % set(gca,'FontSize', 20);
% % file4 =  char(strcat(path,'\DB_dom_4'));
% % saveas(f10,file4,'png')
% % saveas(f10,file4,'fig')

% yamaguchi components

folderName = strcat(FolderName,'\','Yamaguchi4_Y4R_Dbl.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
Y4R_dbl = fread(fileID,[a b],'float32');
Y4R_dbl = Y4R_dbl';

folderName = strcat(FolderName,'\','Yamaguchi4_Y4R_Odd.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
Y4R_odd = fread(fileID,[a b],'float32');
Y4R_odd = Y4R_odd';

folderName = strcat(FolderName,'\','Yamaguchi4_Y4R_Vol.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
Y4R_vol = fread(fileID,[a b],'float32');
Y4R_vol = Y4R_vol';

folderName = strcat(FolderName,'\','Yamaguchi4_Y4R_Hlx.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
Y4R_hlx = fread(fileID,[a b],'float32');
Y4R_hlx = Y4R_hlx';


% dom_class_Y4R_1 = zeros(nrow,ncol);
% dom_class_Y4R_2 = zeros(nrow,ncol);
% dom_class_Y4R_3 = zeros(nrow,ncol);
% dom_class_Y4R_4 = zeros(nrow,ncol);
% 
% 
% for ii = 1:nrow
%     for jj = 1:ncol
%         vec = [Y4R_dbl(ii,jj),Y4R_odd(ii,jj),Y4R_vol(ii,jj),Y4R_hlx(ii,jj)];
%         [B,I] = sort(vec,'descend');
% %         dom_class_Y4R_1(ii,jj) = I(1,1);
% %         dom_class_Y4R_2(ii,jj) = I(1,2);
% %         dom_class_Y4R_3(ii,jj) = I(1,3);
% %         dom_class_Y4R_4(ii,jj) = I(1,4);
%         dom = I(1,1);
%         dom_norm = B(1,1)/sum(B);
%         if dom == 1
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_Y4R_1(ii,jj) = 1;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_Y4R_1(ii,jj) = 2;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_Y4R_1(ii,jj) = 3;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_Y4R_1(ii,jj) = 4;
%             end
%         end
%         
%         if dom == 2
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_Y4R_1(ii,jj) = 5;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_Y4R_1(ii,jj) = 6;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_Y4R_1(ii,jj) = 7;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_Y4R_1(ii,jj) = 8;
%             end
%         end
%         
%         if dom == 3
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_Y4R_1(ii,jj) = 9;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_Y4R_1(ii,jj) = 10;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_Y4R_1(ii,jj) = 11;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_Y4R_1(ii,jj) = 12;
%             end
%         end
%         
%         if dom == 4
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_Y4R_1(ii,jj) = 13;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_Y4R_1(ii,jj) = 14;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_Y4R_1(ii,jj) = 15;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_Y4R_1(ii,jj) = 16;
%             end
%         end
%     end
% end
% 
% f7 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% set(gca,'FontSize',20)
% imagesc(dom_class_Y4R_1)
% axis('image');
% axis off;
% 
% mymap = [0.4,0,0;
%     1,0,0;
%     0.88,0.24,0.16;
%     0.80,0.36,0.36;%red
%     0.53,0.1,0.98;
%     0.12,0.56,1;
%     0.48,0.41,0.93;
%     0,0,1;%blue
%     0.68,1,0.18;
%     0.56,0.74,0.56;
%     0.13,0.54,0.13;
%     0,1,0;%green
%     1,0.08,0.57;
%     1,0.41,0.71;
%     1,0.71,0.76;
%     1,0.74,0.79;%pink
%    ];
% 
% colormap(mymap);
% 
% numcolors = 16;
% caxis([1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% cbarHandle = colorbar('YTick',...
% [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% 'YTickLabel',{'Z1','Z2','Z3','Z4','Z5','Z6','Z7','Z8','Z9','Z10','Z11','Z12','Z13','Z14','Z15','Z16'}, 'YLim', [1 numcolors]);
% set(gca,'FontSize', 20);
% file1 =  char(strcat(path,'\Y4R_dom_1'));
% saveas(f7,file1,'png')
% saveas(f7,file1,'fig')
% 
% % f8 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% % set(gca,'FontSize',20)
% % imagesc(dom_class_Y4R_2)
% % axis('image');
% % axis off;
% % 
% % mymap = [1,0,0;
% %     0,0,1;
% %     0,1,0;
% %     1,0.76,0.8;
% %    ];
% % 
% % colormap(mymap);
% % 
% % numcolors = 4;
% % caxis([1 numcolors]);
% % % cbarHandle = colorbar('YTick',...
% % % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',{'Even','Odd','Diff','Hlx'}, 'YLim', [1 numcolors]);
% % set(gca,'FontSize', 20);
% % file2 =  char(strcat(path,'\Y4R_dom_2'));
% % saveas(f8,file2,'png')
% % saveas(f8,file2,'fig')
% % 
% % f9 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% % set(gca,'FontSize',20)
% % imagesc(dom_class_Y4R_3)
% % axis('image');
% % axis off;
% % 
% % mymap = [1,0,0;
% %     0,0,1;
% %     0,1,0;
% %     1,0.76,0.8;
% %    ];
% % 
% % colormap(mymap);
% % 
% % numcolors = 4;
% % caxis([1 numcolors]);
% % % cbarHandle = colorbar('YTick',...
% % % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',{'Even','Odd','Diff','Hlx'}, 'YLim', [1 numcolors]);
% % set(gca,'FontSize', 20);
% % file3 =  char(strcat(path,'\Y4R_dom_3'));
% % saveas(f9,file3,'png')
% % saveas(f9,file3,'fig')
% % 
% % f10 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% % set(gca,'FontSize',20)
% % imagesc(dom_class_Y4R_4)
% % axis('image');
% % axis off;
% % 
% % mymap = [1,0,0;
% %     0,0,1;
% %     0,1,0;
% %     1,0.76,0.8;
% %    ];
% % 
% % colormap(mymap);
% % 
% % numcolors = 4;
% % caxis([1 numcolors]);
% % % cbarHandle = colorbar('YTick',...
% % % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',{'Even','Odd','Diff','Hlx'}, 'YLim', [1 numcolors]);
% % set(gca,'FontSize', 20);
% % file4 =  char(strcat(path,'\Y4R_dom_4'));
% % saveas(f10,file4,'png')
% % saveas(f10,file4,'fig')


im2 = cat(4,Y4R_dbl',Y4R_vol',Y4R_odd',Y4R_hlx');
f_name_allpowers = strcat(['Y4R_Pd_Pv_Ps_Pc','.bin']);
fileandpath_allpowers = strcat([path f_name_allpowers]);
fid_10 = fopen(fileandpath_allpowers,'wb');
fwrite(fid_10,im2, 'float32');
fclose('all');

Y4R_dbl_sub = Y4R_dbl(top_x:bottom_x,top_y:bottom_y);
Y4R_odd_sub = Y4R_odd(top_x:bottom_x,top_y:bottom_y);
Y4R_hlx_sub = Y4R_hlx(top_x:bottom_x,top_y:bottom_y);
Y4R_vol_sub = Y4R_vol(top_x:bottom_x,top_y:bottom_y);

% average
Y4R_dbl_mean = mean(Y4R_dbl_sub(:));
Y4R_odd_mean = mean(Y4R_odd_sub(:));
Y4R_hlx_mean = mean(Y4R_hlx_sub(:));
Y4R_vol_mean = mean(Y4R_vol_sub(:));

% percentages
Y4R_dbl_percent = Y4R_dbl_mean./(Y4R_dbl_mean + Y4R_odd_mean + Y4R_hlx_mean + Y4R_vol_mean);
Y4R_odd_percent = Y4R_odd_mean./(Y4R_dbl_mean + Y4R_odd_mean + Y4R_hlx_mean + Y4R_vol_mean);
Y4R_hlx_percent = Y4R_hlx_mean./(Y4R_dbl_mean + Y4R_odd_mean + Y4R_hlx_mean + Y4R_vol_mean);
Y4R_vol_percent = Y4R_vol_mean./(Y4R_dbl_mean + Y4R_odd_mean + Y4R_hlx_mean + Y4R_vol_mean);

X = [Y4R_dbl_percent,Y4R_odd_percent,Y4R_hlx_percent,Y4R_vol_percent];

fig3 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
labels = {' ',' ',' ',' '};
pie(X,labels);
colormap([1 0 0;      %// red
          0 0 1;      %// blue
          1 0.08 0.6;      %// red_deriv
          0 1 0])  %// green

% user_entry = input(' give file name: ', 's');
user_entry = strcat('Y4R_',user_entry);
file1 =  char(user_entry);
saveas(fig3,file1,'png')
saveas(fig3,file1,'fig')

save_txt = strcat(user_entry,'.txt');
fid = fopen(save_txt,'wt');

for jj = 1:4
    if jj == 1
        str_p = 'y4r_pd';
    end
    if jj == 2
        str_p = 'y4r_ps';
    end
    if jj == 3
        str_p = 'y4r_pc';
    end
    if jj == 4
        str_p = 'y4r_pv';
    end
    fprintf(fid,'%s\t',str_p);
end

fprintf(fid,'\n');

for ii = 1:size(X,1)
    fprintf(fid,'%g\t',X(ii,:).*100);
    fprintf(fid,'\n');
end

fclose(fid);
close all;

%% G4U

folderName = strcat(FolderName,'\','Singh4_G4U1_Dbl.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
Y4R_dbl = fread(fileID,[a b],'float32');
Y4R_dbl = Y4R_dbl';

folderName = strcat(FolderName,'\','Singh4_G4U1_Odd.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
Y4R_odd = fread(fileID,[a b],'float32');
Y4R_odd = Y4R_odd';

folderName = strcat(FolderName,'\','Singh4_G4U1_Vol.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
Y4R_vol = fread(fileID,[a b],'float32');
Y4R_vol = Y4R_vol';

folderName = strcat(FolderName,'\','Singh4_G4U1_Hlx.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
Y4R_hlx = fread(fileID,[a b],'float32');
Y4R_hlx = Y4R_hlx';


%%
% dom_class_G4U_1 = zeros(nrow,ncol);
% dom_class_G4U_2 = zeros(nrow,ncol);
% dom_class_G4U_3 = zeros(nrow,ncol);
% dom_class_G4U_4 = zeros(nrow,ncol);
% 
% 
% for ii = 1:nrow
%     for jj = 1:ncol
%         vec = [Y4R_dbl(ii,jj),Y4R_odd(ii,jj),Y4R_vol(ii,jj),Y4R_hlx(ii,jj)];
%         [B,I] = sort(vec,'descend');
% %         dom_class_G4U_1(ii,jj) = I(1,1);
% %         dom_class_G4U_2(ii,jj) = I(1,2);
% %         dom_class_G4U_3(ii,jj) = I(1,3);
% %         dom_class_G4U_4(ii,jj) = I(1,4);
%         dom = I(1,1);
%         dom_norm = B(1,1)/sum(B);
%         if dom == 1
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_G4U_1(ii,jj) = 1;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_G4U_1(ii,jj) = 2;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_G4U_1(ii,jj) = 3;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_G4U_1(ii,jj) = 4;
%             end
%         end
%         
%         if dom == 2
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_G4U_1(ii,jj) = 5;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_G4U_1(ii,jj) = 6;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_G4U_1(ii,jj) = 7;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_G4U_1(ii,jj) = 8;
%             end
%         end
%         
%         if dom == 3
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_G4U_1(ii,jj) = 9;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_G4U_1(ii,jj) = 10;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_G4U_1(ii,jj) = 11;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_G4U_1(ii,jj) = 12;
%             end
%         end
%         
%         if dom == 4
%             if dom_norm >= 0 && dom_norm < 0.25
%                 dom_class_G4U_1(ii,jj) = 13;
%             end
%             if dom_norm >= 0.25 && dom_norm < 0.5
%                 dom_class_G4U_1(ii,jj) = 14;
%             end
%             if dom_norm >= 0.5 && dom_norm < 0.75
%                 dom_class_G4U_1(ii,jj) = 15;
%             end
%             if dom_norm >= 0.75 && dom_norm <= 1.0
%                 dom_class_G4U_1(ii,jj) = 16;
%             end
%         end
%     end
% end
% 
% f7 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% set(gca,'FontSize',20)
% imagesc(dom_class_G4U_1)
% axis('image');
% axis off;
% 
% mymap = [0.4,0,0;
%     1,0,0;
%     0.88,0.24,0.16;
%     0.80,0.36,0.36;%red
%     0.53,0.1,0.98;
%     0.12,0.56,1;
%     0.48,0.41,0.93;
%     0,0,1;%blue
%     0.68,1,0.18;
%     0.56,0.74,0.56;
%     0.13,0.54,0.13;
%     0,1,0;%green
%     1,0.08,0.57;
%     1,0.41,0.71;
%     1,0.71,0.76;
%     1,0.74,0.79;%pink
%    ];
% 
% colormap(mymap);
% 
% numcolors = 16;
% caxis([1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% cbarHandle = colorbar('YTick',...
% [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% 'YTickLabel',{'Z1','Z2','Z3','Z4','Z5','Z6','Z7','Z8','Z9','Z10','Z11','Z12','Z13','Z14','Z15','Z16'}, 'YLim', [1 numcolors]);
% set(gca,'FontSize', 20);
% file1 =  char(strcat(path,'\G4U_dom_1'));
% saveas(f7,file1,'png')
% saveas(f7,file1,'fig')
% 
% % f8 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% % set(gca,'FontSize',20)
% % imagesc(dom_class_G4U_2)
% % axis('image');
% % axis off;
% % 
% % mymap = [1,0,0;
% %     0,0,1;
% %     0,1,0;
% %     1,0.76,0.8;
% %    ];
% % 
% % colormap(mymap);
% % 
% % numcolors = 4;
% % caxis([1 numcolors]);
% % % cbarHandle = colorbar('YTick',...
% % % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',{'Even','Odd','Diff','Hlx'}, 'YLim', [1 numcolors]);
% % set(gca,'FontSize', 20);
% % file2 =  char(strcat(path,'\G4U_dom_2'));
% % saveas(f8,file2,'png')
% % saveas(f8,file2,'fig')
% % 
% % f9 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% % set(gca,'FontSize',20)
% % imagesc(dom_class_G4U_3)
% % axis('image');
% % axis off;
% % 
% % mymap = [1,0,0;
% %     0,0,1;
% %     0,1,0;
% %     1,0.76,0.8;
% %    ];
% % 
% % colormap(mymap);
% % 
% % numcolors = 4;
% % caxis([1 numcolors]);
% % % cbarHandle = colorbar('YTick',...
% % % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',{'Even','Odd','Diff','Hlx'}, 'YLim', [1 numcolors]);
% % set(gca,'FontSize', 20);
% % file3 =  char(strcat(path,'\G4U_dom_3'));
% % saveas(f9,file3,'png')
% % saveas(f9,file3,'fig')
% % 
% % f10 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
% % set(gca,'FontSize',20)
% % imagesc(dom_class_G4U_4)
% % axis('image');
% % axis off;
% % 
% % mymap = [1,0,0;
% %     0,0,1;
% %     0,1,0;
% %     1,0.76,0.8;
% %    ];
% % 
% % colormap(mymap);
% % 
% % numcolors = 4;
% % caxis([1 numcolors]);
% % % cbarHandle = colorbar('YTick',...
% % % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % % 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
% % cbarHandle = colorbar('YTick',...
% % [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% % 'YTickLabel',{'Even','Odd','Diff','Hlx'}, 'YLim', [1 numcolors]);
% % set(gca,'FontSize', 20);
% % file4 =  char(strcat(path,'\G4U_dom_4'));
% % saveas(f10,file4,'png')
% % saveas(f10,file4,'fig')
%%

im3 = cat(4,Y4R_dbl',Y4R_vol',Y4R_odd',Y4R_hlx');
f_name_allpowers = strcat(['G4U_Pd_Pv_Ps_Pc','.bin']);
fileandpath_allpowers = strcat([path f_name_allpowers]);
fid_11 = fopen(fileandpath_allpowers,'wb');
fwrite(fid_11,im3, 'float32');
fclose('all');

Y4R_dbl_sub = Y4R_dbl(top_x:bottom_x,top_y:bottom_y);
Y4R_odd_sub = Y4R_odd(top_x:bottom_x,top_y:bottom_y);
Y4R_hlx_sub = Y4R_hlx(top_x:bottom_x,top_y:bottom_y);
Y4R_vol_sub = Y4R_vol(top_x:bottom_x,top_y:bottom_y);

% average
Y4R_dbl_mean = mean(Y4R_dbl_sub(:));
Y4R_odd_mean = mean(Y4R_odd_sub(:));
Y4R_hlx_mean = mean(Y4R_hlx_sub(:));
Y4R_vol_mean = mean(Y4R_vol_sub(:));

% percentages
Y4R_dbl_percent = Y4R_dbl_mean./(Y4R_dbl_mean + Y4R_odd_mean + Y4R_hlx_mean + Y4R_vol_mean);
Y4R_odd_percent = Y4R_odd_mean./(Y4R_dbl_mean + Y4R_odd_mean + Y4R_hlx_mean + Y4R_vol_mean);
Y4R_hlx_percent = Y4R_hlx_mean./(Y4R_dbl_mean + Y4R_odd_mean + Y4R_hlx_mean + Y4R_vol_mean);
Y4R_vol_percent = Y4R_vol_mean./(Y4R_dbl_mean + Y4R_odd_mean + Y4R_hlx_mean + Y4R_vol_mean);

X = [Y4R_dbl_percent,Y4R_odd_percent,Y4R_hlx_percent,Y4R_vol_percent];

fig3 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
labels = {' ',' ',' ',' '};
pie(X,labels);
colormap([1 0 0;      %// red
          0 0 1;      %// blue
          1 0.08 0.6;      %// red_deriv
          0 1 0])  %// green

% user_entry = input(' give file name: ', 's');
user_entry = strcat('G4U_',user_entry);
file1 =  char(user_entry);
saveas(fig3,file1,'png')
saveas(fig3,file1,'fig')

save_txt = strcat(user_entry,'.txt');
fid = fopen(save_txt,'wt');

for jj = 1:4
    if jj == 1
        str_p = 'g4u_pd';
    end
    if jj == 2
        str_p = 'g4u_ps';
    end
    if jj == 3
        str_p = 'g4u_pc';
    end
    if jj == 4
        str_p = 'g4u_pv';
    end
    fprintf(fid,'%s\t',str_p);
end

fprintf(fid,'\n');

for ii = 1:size(X,1)
    fprintf(fid,'%g\t',X(ii,:).*100);
    fprintf(fid,'\n');
end

fclose(fid);
close all;
%%

warning('off');
clc;
close all;
fontSize = 16;
imshow(theta_val_f, []);
axis on;
title('Original Grayscale Image', 'FontSize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

% message = sprintf('Left click to select top pixel');
% uiwait(msgbox(message));
hFH = imfreehand();
% Create a binary image ("mask") from the ROI object.
% binaryImage = hFH.createMask();
xy = hFH.getPosition;
top_pos = xy(1,:);
% message = sprintf('Left click to select bottom pixel');
% uiwait(msgbox(message));
% hFH = imfreehand();
% Create a binary image ("mask") from the ROI object.
% binaryImage1 = hFH.createMask();
% xy = hFH.getPosition;
% bottom_pos = xy(1,:);
% ROI x and y co-ordinates
top_x = top_pos(1,1);
top_y = top_pos(1,2);
top_pos = int32([top_x top_y]);
top_y = top_pos(1,1);
top_x = top_pos(1,2);
% h = drawcircle('Center',[top_x,top_y],'Radius', 6,'StripeColor','red');

% bottom_x = bottom_pos(1,1);
% bottom_y = bottom_pos(1,2);
% bottom_pos = int32([bottom_x bottom_y]);
% bottom_y = bottom_pos(1,1);
% bottom_x = bottom_pos(1,2);

% h1 = drawcircle('Center',[bottom_x,bottom_y],'Radius', 6,'StripeColor','red');
% subsampled as per ROI
t11r = T11(top_x,top_y);
t12r = T12(top_x,top_y);
t13r = T13(top_x,top_y);
t21r = T21(top_x,top_y);
t22r = T22(top_x,top_y);
t23r = T23(top_x,top_y);
t31r = T31(top_x,top_y);
t32r = T32(top_x,top_y);
t33r = T33(top_x,top_y);

T_r = [t11r, t12r, t13r; t21r, t22r, t23r; t31r, t32r, t33r];
theta_r = theta_val_f(top_x,top_y);
tau_r = tau_val_f(top_x,top_y);

disp('T:')
disp(T_r)
disp('theta:')
disp(theta_r)
disp('tau:')
disp(tau_r)
close all;