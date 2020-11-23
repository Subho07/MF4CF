[filename, path] = uigetfile('*.txt', 'Select Full Pol config file');
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
theta_dbl_f = zeros(nrow, ncol); % theta_dbl


disp('load complete');

%%
wsi = 7;
disp('Taking window size as 7')
wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column

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
        
        

        s0_f = trace(T);

        
        val1 = (4*dop_f*k11_f*k44_f)./(k44_f^2 - (1 + 4*dop_f^2)*k11_f^2);
        val2 = abs(k14_f)./(k11_f);
       
        theta_f = atand(val1);
        tau_f = atand(val2);
        theta_val_f(ii,jj) = theta_f;
        tau_val_f(ii,jj) = tau_f;
         
        pc_f(ii,jj) = dop_f.*s0_f.*sind(2.*tau_f);
        pv_f(ii,jj) = (1-dop_f).*s0_f;
        res_pow = s0_f - (pc_f(ii,jj) + pv_f(ii,jj));
        ps_f(ii,jj) = (res_pow/2).*(1+sind(2*theta_f));
        pd_f(ii,jj) = (res_pow/2).*(1-sind(2*theta_f));
  end
  fprintf('row: %d\n',ii);
end
%%
path = FolderName;
if ~exist(strcat(path,'MF4C'), 'dir')
    mkdir(strcat(path,'MF4C'));
end
fclose('all');
FolderName = path;
Fold = strcat(FolderName,'MF4C\');

path = Fold;
cd(path);
%%

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

fName = 'theta_dbl';
f_name_thetaDBL = strcat('\', char(fName), '.bin');
fileandpath_thetaDBL=strcat([path f_name_thetaDBL]);
fid_01 = fopen(fileandpath_thetaDBL,'wb');
fwrite(fid_01,theta_dbl_f', 'float32');
hdrwrite_envi(fName, path, nrow, ncol);
disp('theta_dbl bin and hdr files written');

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

fclose('all');
