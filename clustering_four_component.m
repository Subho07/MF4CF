clc;
clear all;
close all;

[filename, path] = uigetfile('*.bin', 'Select pd');
cd(path)
FolderName = path;

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

folderName = strcat(FolderName,'\',filename);
disp(folderName);
fileID = fopen(folderName,'rb');
pd = fread(fileID,[a b],'float32');
pd = pd';

[filename, path] = uigetfile('*.bin', 'Select ps');

folderName = strcat(FolderName,'\',filename);
disp(folderName);
fileID = fopen(folderName,'rb');
ps = fread(fileID,[a b],'float32');
ps = ps';

[filename, path] = uigetfile('*.bin', 'Select pv');

folderName = strcat(FolderName,'\',filename);
disp(folderName);
fileID = fopen(folderName,'rb');
pv = fread(fileID,[a b],'float32');
pv = pv';

[filename, path] = uigetfile('*.bin', 'Select pc');

folderName = strcat(FolderName,'\',filename);
disp(folderName);
fileID = fopen(folderName,'rb');
pc = fread(fileID,[a b],'float32');
pc = pc';
%%
dom_class_DB_1 = zeros(nrow,ncol);

frac_val = 0.5;

for ii = 1:nrow
    for jj = 1:ncol
        vec = [pd(ii,jj),ps(ii,jj),pv(ii,jj),pc(ii,jj)];
        [B,I] = sort(vec,'descend');

        dom1 = I(1,1);
        dom2 = I(1,2);
        dom3 = I(1,3);
        dom4 = I(1,4);
        
        dom_pow = B(1,1);
        tot_pow = sum(B);
        frac_dom = dom_pow/tot_pow;
        
        if frac_dom <= frac_val
            dom_class_DB_1(ii,jj) = 99;
        %
        elseif (dom1 == 1) && (dom2 == 2) && (dom3 == 3) && (dom4 ==4)
            dom_class_DB_1(ii,jj) = 1;
        elseif (dom1 == 1) && (dom2 == 2) && (dom3 == 4) && (dom4 ==3)
            dom_class_DB_1(ii,jj) = 2;
        elseif (dom1 == 1) && (dom2 == 3) && (dom3 == 2) && (dom4 ==4)
            dom_class_DB_1(ii,jj) = 3;
        elseif (dom1 == 1) && (dom2 == 3) && (dom3 == 4) && (dom4 ==2)
            dom_class_DB_1(ii,jj) = 4;
        elseif (dom1 == 1) && (dom2 == 4) && (dom3 == 2) && (dom4 ==3)
            dom_class_DB_1(ii,jj) = 5;
        elseif (dom1 == 1) && (dom2 == 4) && (dom3 == 3) && (dom4 ==2)
            dom_class_DB_1(ii,jj) = 6;
        %
        elseif (dom1 == 2) && (dom2 == 1) && (dom3 == 3) && (dom4 ==4)
            dom_class_DB_1(ii,jj) = 7;
        elseif (dom1 == 2) && (dom2 == 1) && (dom3 == 4) && (dom4 ==3)
            dom_class_DB_1(ii,jj) = 8;
        elseif (dom1 == 2) && (dom2 == 3) && (dom3 == 1) && (dom4 ==4)
            dom_class_DB_1(ii,jj) = 9;
        elseif (dom1 == 2) && (dom2 == 3) && (dom3 == 4) && (dom4 ==1)
            dom_class_DB_1(ii,jj) = 10;
        elseif (dom1 == 2) && (dom2 == 4) && (dom3 == 1) && (dom4 ==3)
            dom_class_DB_1(ii,jj) = 11;
        elseif (dom1 == 2) && (dom2 == 4) && (dom3 == 3) && (dom4 ==1)
            dom_class_DB_1(ii,jj) = 12;
        %
        elseif (dom1 == 3) && (dom2 == 1) && (dom3 == 2) && (dom4 ==4)
            dom_class_DB_1(ii,jj) = 13;
        elseif (dom1 == 3) && (dom2 == 1) && (dom3 == 4) && (dom4 ==2)
            dom_class_DB_1(ii,jj) = 14;
        elseif (dom1 == 3) && (dom2 == 2) && (dom3 == 1) && (dom4 ==4)
            dom_class_DB_1(ii,jj) = 15;
        elseif (dom1 == 3) && (dom2 == 2) && (dom3 == 4) && (dom4 ==1)
            dom_class_DB_1(ii,jj) = 16;
        elseif (dom1 == 3) && (dom2 == 4) && (dom3 == 1) && (dom4 ==2)
            dom_class_DB_1(ii,jj) = 17;
        elseif (dom1 == 3) && (dom2 == 4) && (dom3 == 2) && (dom4 ==1)
            dom_class_DB_1(ii,jj) = 18;
        %
        elseif (dom1 == 4) && (dom2 == 1) && (dom3 == 2) && (dom4 ==3)
            dom_class_DB_1(ii,jj) = 19;
        elseif (dom1 == 4) && (dom2 == 1) && (dom3 == 3) && (dom4 ==2)
            dom_class_DB_1(ii,jj) = 20;
        elseif (dom1 == 4) && (dom2 == 2) && (dom3 == 1) && (dom4 ==3)
            dom_class_DB_1(ii,jj) = 21;
        elseif (dom1 == 4) && (dom2 == 2) && (dom3 == 3) && (dom4 ==1)
            dom_class_DB_1(ii,jj) = 22;
        elseif (dom1 == 4) && (dom2 == 3) && (dom3 == 1) && (dom4 ==2)
            dom_class_DB_1(ii,jj) = 23;
        elseif (dom1 == 4) && (dom2 == 3) && (dom3 == 2) && (dom4 ==1)
            dom_class_DB_1(ii,jj) = 24;
        end
            
    end
    fprintf('row: %d\n',ii);
end

numel(find(dom_class_DB_1==99))
pause(2);

mean_pd_1 = sum(pd(dom_class_DB_1==1))/numel(pd(dom_class_DB_1==1));
mean_pd_2 = sum(pd(dom_class_DB_1==2))/numel(pd(dom_class_DB_1==2));
mean_pd_3 = sum(pd(dom_class_DB_1==3))/numel(pd(dom_class_DB_1==3));
mean_pd_4 = sum(pd(dom_class_DB_1==4))/numel(pd(dom_class_DB_1==4));
mean_pd_5 = sum(pd(dom_class_DB_1==5))/numel(pd(dom_class_DB_1==5));
mean_pd_6 = sum(pd(dom_class_DB_1==6))/numel(pd(dom_class_DB_1==6));
mean_pd_7 = sum(pd(dom_class_DB_1==7))/numel(pd(dom_class_DB_1==7));
mean_pd_8 = sum(pd(dom_class_DB_1==8))/numel(pd(dom_class_DB_1==8));
mean_pd_9 = sum(pd(dom_class_DB_1==9))/numel(pd(dom_class_DB_1==9));
mean_pd_10 = sum(pd(dom_class_DB_1==10))/numel(pd(dom_class_DB_1==10));
mean_pd_11 = sum(pd(dom_class_DB_1==11))/numel(pd(dom_class_DB_1==11));
mean_pd_12 = sum(pd(dom_class_DB_1==12))/numel(pd(dom_class_DB_1==12));
mean_pd_13 = sum(pd(dom_class_DB_1==13))/numel(pd(dom_class_DB_1==13));
mean_pd_14 = sum(pd(dom_class_DB_1==14))/numel(pd(dom_class_DB_1==14));
mean_pd_15 = sum(pd(dom_class_DB_1==15))/numel(pd(dom_class_DB_1==15));
mean_pd_16 = sum(pd(dom_class_DB_1==16))/numel(pd(dom_class_DB_1==16));
mean_pd_17 = sum(pd(dom_class_DB_1==17))/numel(pd(dom_class_DB_1==17));
mean_pd_18 = sum(pd(dom_class_DB_1==18))/numel(pd(dom_class_DB_1==18));
mean_pd_19 = sum(pd(dom_class_DB_1==19))/numel(pd(dom_class_DB_1==19));
mean_pd_20 = sum(pd(dom_class_DB_1==20))/numel(pd(dom_class_DB_1==20));
mean_pd_21 = sum(pd(dom_class_DB_1==21))/numel(pd(dom_class_DB_1==21));
mean_pd_22 = sum(pd(dom_class_DB_1==22))/numel(pd(dom_class_DB_1==22));
mean_pd_23 = sum(pd(dom_class_DB_1==23))/numel(pd(dom_class_DB_1==23));
mean_pd_24 = sum(pd(dom_class_DB_1==24))/numel(pd(dom_class_DB_1==24));
%
mean_ps_1 = sum(ps(dom_class_DB_1==1))/numel(ps(dom_class_DB_1==1));
mean_ps_2 = sum(ps(dom_class_DB_1==2))/numel(ps(dom_class_DB_1==2));
mean_ps_3 = sum(ps(dom_class_DB_1==3))/numel(ps(dom_class_DB_1==3));
mean_ps_4 = sum(ps(dom_class_DB_1==4))/numel(ps(dom_class_DB_1==4));
mean_ps_5 = sum(ps(dom_class_DB_1==5))/numel(ps(dom_class_DB_1==5));
mean_ps_6 = sum(ps(dom_class_DB_1==6))/numel(ps(dom_class_DB_1==6));
mean_ps_7 = sum(ps(dom_class_DB_1==7))/numel(ps(dom_class_DB_1==7));
mean_ps_8 = sum(ps(dom_class_DB_1==8))/numel(ps(dom_class_DB_1==8));
mean_ps_9 = sum(ps(dom_class_DB_1==9))/numel(ps(dom_class_DB_1==9));
mean_ps_10 = sum(ps(dom_class_DB_1==10))/numel(ps(dom_class_DB_1==10));
mean_ps_11 = sum(ps(dom_class_DB_1==11))/numel(ps(dom_class_DB_1==11));
mean_ps_12 = sum(ps(dom_class_DB_1==12))/numel(ps(dom_class_DB_1==12));
mean_ps_13 = sum(ps(dom_class_DB_1==13))/numel(ps(dom_class_DB_1==13));
mean_ps_14 = sum(ps(dom_class_DB_1==14))/numel(ps(dom_class_DB_1==14));
mean_ps_15 = sum(ps(dom_class_DB_1==15))/numel(ps(dom_class_DB_1==15));
mean_ps_16 = sum(ps(dom_class_DB_1==16))/numel(ps(dom_class_DB_1==16));
mean_ps_17 = sum(ps(dom_class_DB_1==17))/numel(ps(dom_class_DB_1==17));
mean_ps_18 = sum(ps(dom_class_DB_1==18))/numel(ps(dom_class_DB_1==18));
mean_ps_19 = sum(ps(dom_class_DB_1==19))/numel(ps(dom_class_DB_1==19));
mean_ps_20 = sum(ps(dom_class_DB_1==20))/numel(ps(dom_class_DB_1==20));
mean_ps_21 = sum(ps(dom_class_DB_1==21))/numel(ps(dom_class_DB_1==21));
mean_ps_22 = sum(ps(dom_class_DB_1==22))/numel(ps(dom_class_DB_1==22));
mean_ps_23 = sum(ps(dom_class_DB_1==23))/numel(ps(dom_class_DB_1==23));
mean_ps_24 = sum(ps(dom_class_DB_1==24))/numel(ps(dom_class_DB_1==24));
%
mean_pv_1 = sum(pv(dom_class_DB_1==1))/numel(pv(dom_class_DB_1==1));
mean_pv_2 = sum(pv(dom_class_DB_1==2))/numel(pv(dom_class_DB_1==2));
mean_pv_3 = sum(pv(dom_class_DB_1==3))/numel(pv(dom_class_DB_1==3));
mean_pv_4 = sum(pv(dom_class_DB_1==4))/numel(pv(dom_class_DB_1==4));
mean_pv_5 = sum(pv(dom_class_DB_1==5))/numel(pv(dom_class_DB_1==5));
mean_pv_6 = sum(pv(dom_class_DB_1==6))/numel(pv(dom_class_DB_1==6));
mean_pv_7 = sum(pv(dom_class_DB_1==7))/numel(pv(dom_class_DB_1==7));
mean_pv_8 = sum(pv(dom_class_DB_1==8))/numel(pv(dom_class_DB_1==8));
mean_pv_9 = sum(pv(dom_class_DB_1==9))/numel(pv(dom_class_DB_1==9));
mean_pv_10 = sum(pv(dom_class_DB_1==10))/numel(pv(dom_class_DB_1==10));
mean_pv_11 = sum(pv(dom_class_DB_1==11))/numel(pv(dom_class_DB_1==11));
mean_pv_12 = sum(pv(dom_class_DB_1==12))/numel(pv(dom_class_DB_1==12));
mean_pv_13 = sum(pv(dom_class_DB_1==13))/numel(pv(dom_class_DB_1==13));
mean_pv_14 = sum(pv(dom_class_DB_1==14))/numel(pv(dom_class_DB_1==14));
mean_pv_15 = sum(pv(dom_class_DB_1==15))/numel(pv(dom_class_DB_1==15));
mean_pv_16 = sum(pv(dom_class_DB_1==16))/numel(pv(dom_class_DB_1==16));
mean_pv_17 = sum(pv(dom_class_DB_1==17))/numel(pv(dom_class_DB_1==17));
mean_pv_18 = sum(pv(dom_class_DB_1==18))/numel(pv(dom_class_DB_1==18));
mean_pv_19 = sum(pv(dom_class_DB_1==19))/numel(pv(dom_class_DB_1==19));
mean_pv_20 = sum(pv(dom_class_DB_1==20))/numel(pv(dom_class_DB_1==20));
mean_pv_21 = sum(pv(dom_class_DB_1==21))/numel(pv(dom_class_DB_1==21));
mean_pv_22 = sum(pv(dom_class_DB_1==22))/numel(pv(dom_class_DB_1==22));
mean_pv_23 = sum(pv(dom_class_DB_1==23))/numel(pv(dom_class_DB_1==23));
mean_pv_24 = sum(pv(dom_class_DB_1==24))/numel(pv(dom_class_DB_1==24));
%
mean_pc_1 = sum(pc(dom_class_DB_1==1))/numel(pc(dom_class_DB_1==1));
mean_pc_2 = sum(pc(dom_class_DB_1==2))/numel(pc(dom_class_DB_1==2));
mean_pc_3 = sum(pc(dom_class_DB_1==3))/numel(pc(dom_class_DB_1==3));
mean_pc_4 = sum(pc(dom_class_DB_1==4))/numel(pc(dom_class_DB_1==4));
mean_pc_5 = sum(pc(dom_class_DB_1==5))/numel(pc(dom_class_DB_1==5));
mean_pc_6 = sum(pc(dom_class_DB_1==6))/numel(pc(dom_class_DB_1==6));
mean_pc_7 = sum(pc(dom_class_DB_1==7))/numel(pc(dom_class_DB_1==7));
mean_pc_8 = sum(pc(dom_class_DB_1==8))/numel(pc(dom_class_DB_1==8));
mean_pc_9 = sum(pc(dom_class_DB_1==9))/numel(pc(dom_class_DB_1==9));
mean_pc_10 = sum(pc(dom_class_DB_1==10))/numel(pc(dom_class_DB_1==10));
mean_pc_11 = sum(pc(dom_class_DB_1==11))/numel(pc(dom_class_DB_1==11));
mean_pc_12 = sum(pc(dom_class_DB_1==12))/numel(pc(dom_class_DB_1==12));
mean_pc_13 = sum(pc(dom_class_DB_1==13))/numel(pc(dom_class_DB_1==13));
mean_pc_14 = sum(pc(dom_class_DB_1==14))/numel(pc(dom_class_DB_1==14));
mean_pc_15 = sum(pc(dom_class_DB_1==15))/numel(pc(dom_class_DB_1==15));
mean_pc_16 = sum(pc(dom_class_DB_1==16))/numel(pc(dom_class_DB_1==16));
mean_pc_17 = sum(pc(dom_class_DB_1==17))/numel(pc(dom_class_DB_1==17));
mean_pc_18 = sum(pc(dom_class_DB_1==18))/numel(pc(dom_class_DB_1==18));
mean_pc_19 = sum(pc(dom_class_DB_1==19))/numel(pc(dom_class_DB_1==19));
mean_pc_20 = sum(pc(dom_class_DB_1==20))/numel(pc(dom_class_DB_1==20));
mean_pc_21 = sum(pc(dom_class_DB_1==21))/numel(pc(dom_class_DB_1==21));
mean_pc_22 = sum(pc(dom_class_DB_1==22))/numel(pc(dom_class_DB_1==22));
mean_pc_23 = sum(pc(dom_class_DB_1==23))/numel(pc(dom_class_DB_1==23));
mean_pc_24 = sum(pc(dom_class_DB_1==24))/numel(pc(dom_class_DB_1==24));
%
vec_1 = [mean_pd_1,mean_ps_1,mean_pv_1,mean_pc_1];
vec_2 = [mean_pd_2,mean_ps_2,mean_pv_2,mean_pc_2];
vec_3 = [mean_pd_3,mean_ps_3,mean_pv_3,mean_pc_3];
vec_4 = [mean_pd_4,mean_ps_4,mean_pv_4,mean_pc_4];
vec_5 = [mean_pd_5,mean_ps_5,mean_pv_5,mean_pc_5];
vec_6 = [mean_pd_6,mean_ps_6,mean_pv_6,mean_pc_6];
vec_7 = [mean_pd_7,mean_ps_7,mean_pv_7,mean_pc_7];
vec_8 = [mean_pd_8,mean_ps_8,mean_pv_8,mean_pc_8];
vec_9 = [mean_pd_9,mean_ps_9,mean_pv_9,mean_pc_9];
vec_10 = [mean_pd_10,mean_ps_10,mean_pv_10,mean_pc_10];
vec_11 = [mean_pd_11,mean_ps_11,mean_pv_11,mean_pc_11];
vec_12 = [mean_pd_12,mean_ps_12,mean_pv_12,mean_pc_12];
vec_13 = [mean_pd_13,mean_ps_13,mean_pv_13,mean_pc_13];
vec_14 = [mean_pd_14,mean_ps_14,mean_pv_14,mean_pc_14];
vec_15 = [mean_pd_15,mean_ps_15,mean_pv_15,mean_pc_15];
vec_16 = [mean_pd_16,mean_ps_16,mean_pv_16,mean_pc_16];
vec_17 = [mean_pd_17,mean_ps_17,mean_pv_17,mean_pc_17];
vec_18 = [mean_pd_18,mean_ps_18,mean_pv_18,mean_pc_18];
vec_19 = [mean_pd_19,mean_ps_19,mean_pv_19,mean_pc_19];
vec_20 = [mean_pd_20,mean_ps_20,mean_pv_20,mean_pc_20];
vec_21 = [mean_pd_21,mean_ps_21,mean_pv_21,mean_pc_21];
vec_22 = [mean_pd_22,mean_ps_22,mean_pv_22,mean_pc_22];
vec_23 = [mean_pd_23,mean_ps_23,mean_pv_23,mean_pc_23];
vec_24 = [mean_pd_24,mean_ps_24,mean_pv_24,mean_pc_24];

for ii = 1:nrow
    for jj = 1:ncol
        if (dom_class_DB_1(ii,jj) == 99)
            pd_mix = pd(ii,jj);
            ps_mix = ps(ii,jj);
            pv_mix = pv(ii,jj);
            pc_mix = pc(ii,jj);
            vec_mix = [pd_mix,ps_mix,pv_mix,pc_mix];
            [B_mix,I_mix] = sort(vec_mix,'descend');
            
            if I_mix(1,1) == 1
                dist_1 = norm(vec_1 - vec_mix);
                dist_2 = norm(vec_2 - vec_mix);
                dist_3 = norm(vec_3 - vec_mix);
                dist_4 = norm(vec_4 - vec_mix);
                dist_5 = norm(vec_5 - vec_mix);
                dist_6 = norm(vec_6 - vec_mix);
                vec1 = [dist_1,dist_2,dist_3,dist_4,dist_5,dist_6];
                index = [1,2,3,4,5,6];
                [~,I1] = sort(vec1,'ascend');
                dom_class_DB_1(ii,jj) = index(1,I1(1,1));
            end
            if I_mix(1,1) == 2
                dist_7 = norm(vec_7 - vec_mix);
                dist_8 = norm(vec_8 - vec_mix);
                dist_9 = norm(vec_9 - vec_mix);
                dist_10 = norm(vec_10 - vec_mix);
                dist_11 = norm(vec_11 - vec_mix);
                dist_12 = norm(vec_12 - vec_mix);
                vec1 = [dist_7,dist_8,dist_9,dist_10,dist_11,dist_12];
                index = [7,8,9,10,11,12];
                [~,I1] = sort(vec1,'ascend');
                dom_class_DB_1(ii,jj) = index(1,I1(1,1));
            end
            if I_mix(1,1) == 3
                dist_13 = norm(vec_13 - vec_mix);
                dist_14 = norm(vec_14 - vec_mix);
                dist_15 = norm(vec_15 - vec_mix);
                dist_16 = norm(vec_16 - vec_mix);
                dist_17 = norm(vec_17 - vec_mix);
                dist_18 = norm(vec_18 - vec_mix);
                vec1 = [dist_13,dist_14,dist_15,dist_16,dist_17,dist_18];
                index = [13,14,15,16,17,18];
                [~,I1] = sort(vec1,'ascend');
                dom_class_DB_1(ii,jj) = index(1,I1(1,1));
            end
            if I_mix(1,1) == 4
                dist_19 = norm(vec_19 - vec_mix);
                dist_20 = norm(vec_20 - vec_mix);
                dist_21 = norm(vec_21 - vec_mix);
                dist_22 = norm(vec_22 - vec_mix);
                dist_23 = norm(vec_23 - vec_mix);
                dist_24 = norm(vec_24 - vec_mix);
                vec1 = [dist_19,dist_20,dist_21,dist_22,dist_23,dist_24];
                index = [19,20,21,22,23,24];
                [~,I1] = sort(vec1,'ascend');
                dom_class_DB_1(ii,jj) = index(1,I1(1,1));
            end
        end
    end
    fprintf('row: %d\n',ii);
end


numel(find(dom_class_DB_1==99))
pause(2);
% imagesc(flipdim(dom_class_DB_1,1));
%%
f7 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
set(gca,'FontSize',20)
imagesc(dom_class_DB_1)
axis('image');
axis off;
% 
% mymap = [0.72,0,0;
%     1,0,0;
%     1,0.22,0.22;
%     1,0.34,0.34;
%     1,0.46,0.46;
%     1,0.64,0.64;%red
%     0,0,0.72;
%     0,0,1;
%     0.22,0.22,1;
%     0.34,0.34,1;
%     0.46,0.46,1;
%     0.64,0.64,1;%blue
%     0,0.72,0;
%     0,1,0;
%     0.22,1,0.22;
%     0.34,1,0.34;
%     0.46,1,0.46;
%     0.64,1,0.64;%green
%     0.77,0.06,0.77;
%     0.84,0.07,0.84;
%     0.93,0.11,0.93;
%     0.94,0.34,0.94;
%     0.96,0.46,0.96;
%     0.97,0.62,0.97;%pink
%    ];
% % 
% colormap(mymap);

hex = ['#dc0032';'#fb2647';'#ff5c63';'#ff8180';'#ffa29e';'#ffc2bd';...
    '#1304d7';'#603adf';'#8961e7';'#ab87ee';'#c9aef4';'#e5d6fa';...
'#00ce00';'#26fa14';'#6bfc52';'#93fe79';'#b2ff9c';'#ceffbd';...
'#da00d8';'#f818fa';'#fd5cfc';'#ff84fd';'#ffa6fe';'#ffc5ff'];
mymap = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
colormap(mymap);

numcolors = 24;
caxis([1 numcolors]);
% cbarHandle = colorbar('YTick',...
% [1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
% 'YTickLabel',int2str([1:numcolors]'), 'YLim', [1 numcolors]);
cbarHandle = colorbar('YTick',...
[1+0.5*(numcolors-1)/numcolors:(numcolors-1)/numcolors:numcolors],...
'YTickLabel',{'Z1','Z2','Z3','Z4','Z5','Z6','Z7','Z8','Z9','Z10','Z11','Z12','Z13','Z14','Z15','Z16'...
'Z17','Z18', 'Z19','Z20','Z21','Z22','Z23','Z24'}, 'YLim', [1 numcolors]);
set(gca,'FontSize', 12);
file1 =  char(strcat(path,'\MF4C_dom_1'));
saveas(f7,file1,'png')
saveas(f7,file1,'fig')

fName = 'MF4C_Clustered';
f_name_pdDP = strcat('\', char(fName), '.bin');
fileandpath_pdDP=strcat([path f_name_pdDP]);
fid_05 = fopen(fileandpath_pdDP,'wb');
fwrite(fid_05,dom_class_DB_1, 'float32');
% hdrwrite_envi(fName, path, nrow, ncol);
disp('cluster file: written');
%%

warning('off');
clc;
close all;
fontSize = 16;
imshow(dom_class_DB_1, []);
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
classes = dom_class_DB_1(top_x:bottom_x,top_y:bottom_y);

classes = classes(:);

[a,b] = hist(classes,unique(classes));

disp(a)
disp(b)
close all;
%%

% warning('off');
% clc;
% close all;
% fontSize = 16;
% imshow(dom_class_DB_1, []);
% axis on;
% title('Original Grayscale Image', 'FontSize', fontSize);
% set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
% 
% % message = sprintf('Left click to select top pixel');
% % uiwait(msgbox(message));
% hFH = imfreehand();
% % Create a binary image ("mask") from the ROI object.
% % binaryImage = hFH.createMask();
% xy = hFH.getPosition;
% top_pos = xy(1,:);
% % ROI x and y co-ordinates
% top_x = top_pos(1,1);
% top_y = top_pos(1,2);
% top_pos = int32([top_x top_y]);
% top_y = top_pos(1,1);
% top_x = top_pos(1,2);
%%
% ii = top_y;
% jj = top_x;
% ii = 623;
% jj = 113;
% 
% vec = [pd(ii,jj),ps(ii,jj),pv(ii,jj),pc(ii,jj)];
% [B,I] = sort(vec,'descend');
% disp(vec);
% disp(B);
% disp(I);
% disp(dom_class_DB_1(ii,jj));
% 
% close all