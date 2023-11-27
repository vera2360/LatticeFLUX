% %% 读取picaso-localization高斯拟合后total3_locs.hdf5数据，进行相位拟合
% % by cxj 
% % 2023.4.21
clc
clear
tic
% %读取子文件,uint16类型，
% DATA_PP=tiffread("C:\Users\15600\Desktop\filament\fila_bk225_1.tif");
% DATA_P0=tiffread("C:\Users\15600\Desktop\filament\fila_bk225_2.tif");
% DATA_PN=tiffread("C:\Users\15600\Desktop\filament\fila_bk225_3.tif");
% DATA_0P=tiffread("C:\Users\15600\Desktop\filament\fila_bk225_4.tif");
% DATA_00=tiffread("C:\Users\15600\Desktop\filament\fila_bk225_5.tif");
% DATA_0N=tiffread("C:\Users\15600\Desktop\filament\fila_bk225_6.tif");
% DATA_NP=tiffread("C:\Users\15600\Desktop\filament\fila_bk225_7.tif");
% DATA_N0=tiffread("C:\Users\15600\Desktop\filament\fila_bk225_8.tif");
% DATA_NN=tiffread("C:\Users\15600\Desktop\filament\fila_bk225_9.tif");
% DATA_total3=tiffread("C:\Users\15600\Desktop\filament\fila_bk225_total.tif");
% toc
load('C:\Users\15600\Desktop\filament\fila_bk500_1.mat');
DATA_PP=fila_lattice_downsample_noise_gain;
load('C:\Users\15600\Desktop\filament\fila_bk500_2.mat');
DATA_P0=fila_lattice_downsample_noise_gain;
load('C:\Users\15600\Desktop\filament\fila_bk500_3.mat');
DATA_PN=fila_lattice_downsample_noise_gain;
load('C:\Users\15600\Desktop\filament\fila_bk500_4.mat');
DATA_0P=fila_lattice_downsample_noise_gain;
load('C:\Users\15600\Desktop\filament\fila_bk500_5.mat');
DATA_00=fila_lattice_downsample_noise_gain;
load('C:\Users\15600\Desktop\filament\fila_bk500_6.mat');
DATA_0N=fila_lattice_downsample_noise_gain;
load('C:\Users\15600\Desktop\filament\fila_bk500_7.mat');
DATA_NP=fila_lattice_downsample_noise_gain;
load('C:\Users\15600\Desktop\filament\fila_bk500_8.mat');
DATA_N0=fila_lattice_downsample_noise_gain;
load('C:\Users\15600\Desktop\filament\fila_bk500_9.mat');
DATA_NN=fila_lattice_downsample_noise_gain;
load('C:\Users\15600\Desktop\filament\fila_lattice_downsample_noise_total.mat');
DATA_TOTAL=fila_lattice_downsample_noise_total;
filename_wformat=...
"C:\Users\15600\Desktop\filament\fila_bk500_total_locs.hdf5";%xy定位文件
info = h5info(filename_wformat);
data = h5read(filename_wformat,'/locs');
x_identify=data.x+ones(size(data.x,1),1);%以pixel为单位,python的矩阵都是从0开始的，所以都要加1
y_identify=data.y+ones(size(data.x,1),1);%以pixel为单位
frame_identify=data.frame+uint32(ones(size(data.frame,1),1));
sx_identify=data.sx;
sy_identify=data.sy;
photons_identify=data.photons;
bg_identify=data.bg;
lpx_identify=data.lpx;
lpy_identify=data.lpy;
ellipticity_identify=data.ellipticity;
net_gradient_identify=data.net_gradient;

%从_total3_locs.yaml中获取locs参数
yaml_file=...
"C:\Users\15600\Desktop\filament\fila_bk500_total_locs.yaml";

yaml_data=yaml.loadFile(yaml_file);
Frame=yaml_data.Frames;
Height=yaml_data.Height;
Width=yaml_data.Width;
boxsize=yaml_data.BoxSize;
pixelsize=yaml_data.Pixelsize;
pitch=330; 
% 
r=floor(boxsize/2);
% r=4;
id1=find(x_identify<3.5 |x_identify>19);id2=find(y_identify<3.5|y_identify>19);
iDATA_D=union(id1,id2);%求两个索引的并集
x_identify(iDATA_D)=[];
y_identify(iDATA_D)=[];
frame_identify(iDATA_D)=[];
sx_identify(iDATA_D)=[];
sy_identify(iDATA_D)=[];
photons_identify(iDATA_D)=[];
bg_identify(iDATA_D)=[];
lpx_identify(iDATA_D)=[];
lpy_identify(iDATA_D)=[];
ellipticity_identify(iDATA_D)=[];
net_gradient_identify(iDATA_D)=[];
% %将所有的locs对应在9个子数据集中的ROI存到4D数组中，第三维是9张，第四维是locs数量
% % % 注意image的坐标系统和矩阵的系统不一致，是互为转置的关系，同时picasso输出的帧数是从0开始计数的，所以需要+1
for i=1:size(frame_identify,1)
    img9(:,:,1,i)=DATA_PP((round(y_identify(i))-r):(round(y_identify(i))+r),(round(x_identify(i))-r):(round(x_identify(i))+r),frame_identify(i));
    img9(:,:,2,i)=DATA_P0((round(y_identify(i))-r):(round(y_identify(i))+r),(round(x_identify(i))-r):(round(x_identify(i))+r),frame_identify(i));
    img9(:,:,3,i)=DATA_PN((round(y_identify(i))-r):(round(y_identify(i))+r),(round(x_identify(i))-r):(round(x_identify(i))+r),frame_identify(i));
    img9(:,:,4,i)=DATA_0P((round(y_identify(i))-r):(round(y_identify(i))+r),(round(x_identify(i))-r):(round(x_identify(i))+r),frame_identify(i));
    img9(:,:,5,i)=DATA_00((round(y_identify(i))-r):(round(y_identify(i))+r),(round(x_identify(i))-r):(round(x_identify(i))+r),frame_identify(i));
    img9(:,:,6,i)=DATA_0N((round(y_identify(i))-r):(round(y_identify(i))+r),(round(x_identify(i))-r):(round(x_identify(i))+r),frame_identify(i));
    img9(:,:,7,i)=DATA_NP((round(y_identify(i))-r):(round(y_identify(i))+r),(round(x_identify(i))-r):(round(x_identify(i))+r),frame_identify(i));
    img9(:,:,8,i)=DATA_N0((round(y_identify(i))-r):(round(y_identify(i))+r),(round(x_identify(i))-r):(round(x_identify(i))+r),frame_identify(i));
    img9(:,:,9,i)=DATA_NN((round(y_identify(i))-r):(round(y_identify(i))+r),(round(x_identify(i))-r):(round(x_identify(i))+r),frame_identify(i));
end
% 
% 执行9个子数据的高斯拟合，得到光子数
tic
for m=1:size(frame_identify,1)
    imgbuf=img9(:,:,:,m);
    cx=x_identify(m);
    cy=y_identify(m);
    sigmax=sx_identify(m);
    sigmay=sy_identify(m);
    fit_9result=fit_9subset(imgbuf,cx,cy,sigmax,sigmay,r,frame_identify(m));    
    % fit_9result=fit_9subset_MLE(imgbuf,cx,cy,sigmax,sigmay,r,frame_identify(m));
    fit(m,:)=fit_9result;
end
toc
    % fit(:,5:13)=fit(:,5:13).*fit(:,3).*fit(:,4).*2.*pi;
    fit(:,1:2)=fit(:,1:2)-ones(size(fit,1),2);
% 计算相位
fit_phase=Cal_phase(fit);
% fit_phase: py_p,py_0,py_n,px_p,px_0,px_n,mdy_p,mdy_0,mdy_n,mdx_p,mdx_0,mdx_n

% % 先调整相位，再计算相位周期
for n=1:size(fit,1)
    fit_phase_adjust(n,1:3)=phase_adjust(fit_phase(n,1:3));
    fit_phase_adjust(n,4:6)=phase_adjust(fit_phase(n,4:6));

    fit_phase_period(n,4:6)=Cal_period(fit(n,1),pitch,pixelsize,fit_phase_adjust(n,4:6));
    fit_phase_period(n,1:3)=Cal_period(fit(n,2),pitch,pixelsize,fit_phase_adjust(n,1:3));%用高斯拟合的位置确定相位周期数

end

fit_sum=[fit(:,1:2)*pixelsize,fit_phase_period/(2*pi)*pitch]; 
% fit_sum: gauss_x,gauss_y,py_p,py_0,py_n,px_p,px_0,px_n

col_Del=any(isnan(fit_sum(:,3:8)),2);
fit_sum_Del=[fit_sum,double(frame_identify),double(photons_identify),double(sx_identify),double(sy_identify),double(bg_identify),...
    double(lpx_identify),double(lpy_identify),double(ellipticity_identify),double(net_gradient_identify)];
%fit_sum_del: x y py1 py2 py3 px1 px2 px3 
fit_sum_Del(col_Del,:)=[];

%% x
roughx=fit_sum_Del(:,1);
phasex1=fit_sum_Del(:,6);
figure(5);
subplot(1,3,1);title("x");hold on
plot(phasex1,roughx,'r.');

label=zeros(size(roughx,1),6);
for u=1:size(roughx,1)
    dist(u,1)=abs(fit_sum_Del(u,7)-fit_sum_Del(u,8));
    dist(u,2)=abs(fit_sum_Del(u,6)-fit_sum_Del(u,8));
    dist(u,3)=abs(fit_sum_Del(u,6)-fit_sum_Del(u,7));
    if max(dist(u,:))-min(dist(u,:))>3
        t_id=find(dist(u,:)==min(dist(u,:)));
        label(u,t_id)=1;
    elseif max(dist(u,:))-min(dist(u,:))<1/2/pi*pitch
        label(u,1:3)=0;
    end
end

iDATA_1=find(label(:,1));
phasex1_fit=phasex1;
phasex1_fit(iDATA_1)=[];
roughx1_fit=roughx;
roughx1_fit(iDATA_1)=[];
px1=polyfit(phasex1_fit,roughx1_fit,1);
refinex1=px1(1).*phasex1+px1(2);
% plot(phasex1_fit,roughx1_fit,'b.');
plot(phasex1,refinex1,'g.');axis square

phasex0=fit_sum_Del(:,7);
subplot(1,3,2);hold on
plot(phasex0,roughx,'r.'); %横坐标是相位，纵坐标是高斯拟合坐标
iDATA_2=find(label(:,2));
phasex0_fit=phasex0;
phasex0_fit(iDATA_2)=[];
roughx0_fit=roughx;
roughx0_fit(iDATA_2)=[];

px0=polyfit(phasex0_fit,roughx0_fit,1); 
refinex0=px0(1).*phasex0+px0(2);
% plot(phasex0_fit,roughx0_fit,'b.');
plot(phasex0,refinex0,'g.');axis square

phasex_1=fit_sum_Del(:,8);
subplot(1,3,3); hold on
plot(phasex_1,roughx,'r.'); %横坐标是相位，纵坐标是高斯拟合坐标
iDATA_3=find(label(:,3));
phasex_1_fit=phasex_1;
phasex_1_fit(iDATA_3)=[];
roughx_1_fit=roughx;
roughx_1_fit(iDATA_3)=[];
px_1=polyfit(phasex_1_fit,roughx_1_fit,1); 
refinex_1=px_1(1).*phasex_1+px_1(2);
% plot(phasex_1_fit,roughx_1_fit,'b.');
plot(phasex_1,refinex_1,'g.');axis square;hold off
%%  y
roughy=fit_sum_Del(:,2);
phasey1=fit_sum_Del(:,3);
figure(7);
subplot(1,3,1);title("y");hold on
plot(phasey1,roughy,'r.'); %横坐标是相位，纵坐标是高斯拟合坐标

for u=1:size(roughy,1)
    dist(u,4)=abs(fit_sum_Del(u,4)-fit_sum_Del(u,5));
    dist(u,5)=abs(fit_sum_Del(u,3)-fit_sum_Del(u,5));
    dist(u,6)=abs(fit_sum_Del(u,3)-fit_sum_Del(u,4));
    if max(dist(u,4:6))-min(dist(u,4:6))>3
        t_id=find(dist(u,4:6)==min(dist(u,4:6)));
        label(u,t_id+3)=1;
    else%if max(dist(u,4:6))-min(dist(u,4:6))<1/2/pi*pitch
        label(u,4:6)=0;
    end
end

iDATA_4=find(label(:,4));
phasey1_fit=phasey1;
phasey1_fit(iDATA_4)=[];
roughy1_fit=roughy;
roughy1_fit(iDATA_4)=[];
% plot(phasey1_fit,roughy1_fit,'b.');
py1=polyfit(phasey1_fit,roughy1_fit,1); 
refiney1=py1(1).*phasey1+py1(2);

plot(phasey1,refiney1,'g.');axis square

phasey0=fit_sum_Del(:,4);
subplot(1,3,2);hold on
plot(phasey0,roughy,'r.'); %横坐标是相位，纵坐标是高斯拟合坐标
iDATA_5=find(label(:,5));
phasey0_fit=phasey0;
phasey0_fit(iDATA_5)=[];
roughy0_fit=roughy;
roughy0_fit(iDATA_5)=[];
% plot(phasey0_fit,roughy0_fit,'b.');
py0=polyfit(phasey0_fit,roughy0_fit,1); 
refiney0=py0(1).*phasey0+py0(2);
plot(phasey0,refiney0,'g.');axis square

phasey_1=fit_sum_Del(:,5);
subplot(1,3,3); hold on
plot(phasey_1,roughy,'r.'); %横坐标是相位，纵坐标是高斯拟合坐标
iDATA_6=find(label(:,6));
phasey_1_fit=phasey_1;
phasey_1_fit(iDATA_6)=[];
roughy_1_fit=roughy;
roughy_1_fit(iDATA_6)=[];
% plot(phasey_1_fit,roughy_1_fit,'b.');
py_1=polyfit(phasey_1_fit,roughy_1_fit,1); 
refiney_1=py_1(1).*phasey_1+py_1(2);
plot(phasey_1,refiney_1,'g.');axis square;hold off

%% tiling，并画晶格调制定位的精度分布
% x相位<=30°，在refinex_1里取，30°<x<=150°，在refinex0里取，150°<x<=270°，在refinex1里取
% y相位同x

refiney_uint=[refiney1,refiney0,refiney_1];
for u=1:size(roughy,1)
    dist_breakpointy(u,1)=abs(refiney_uint(u,2)-refiney_uint(u,3));
    dist_breakpointy(u,2)=abs(refiney_uint(u,1)-refiney_uint(u,3));
    dist_breakpointy(u,3)=abs(refiney_uint(u,1)-refiney_uint(u,2));
    if max(dist_breakpointy(u,1:3))-min(dist_breakpointy(u,1:3))>5 
        breakpoint_id=find(dist_breakpointy(u,1:3)==min(dist_breakpointy(u,1:3))); 
        refiney_uint(u,breakpoint_id)=0;
        refiney_breakpoint(u)=sum(refiney_uint(u,:))/2;
    elseif max(dist_breakpointy(u,1:3))-min(dist_breakpointy(u,1:3))<5
        refiney_breakpoint(u)=sum(refiney_uint(u,:))/3;
    end
end
refinex_uint=[refinex1,refinex0,refinex_1];
for u=1:size(roughy,1)
    dist_breakpoint(u,1)=abs(refinex_uint(u,2)-refinex_uint(u,3));
    dist_breakpoint(u,2)=abs(refinex_uint(u,1)-refinex_uint(u,3));
    dist_breakpoint(u,3)=abs(refinex_uint(u,1)-refinex_uint(u,2));
    if max(dist_breakpoint(u,1:3))-min(dist_breakpoint(u,1:3))>5
        breakpoint_id=find(dist_breakpoint(u,1:3)==min(dist_breakpoint(u,1:3)));
        refinex_uint(u,breakpoint_id)=0;
        refinex_breakpoint(u)=sum(refinex_uint(u,:))/2;
    elseif max(dist_breakpoint(u,1:3))-min(dist_breakpoint(u,1:3))<5
        refinex_breakpoint(u)=sum(refinex_uint(u,:))/3;
    end
end
%%
for i=1:size(refinex0,1)

    if mod(refinex_breakpoint(i),pitch)>=30/360*pitch & mod(refinex_breakpoint(i),pitch)<150/360*pitch %(18.33,91.667)
        refiney_tiling(i,1)=refiney0(i);
    elseif mod(refinex_breakpoint(i),pitch)>=150/360*pitch & mod(refinex_breakpoint(i),pitch)<270/360*pitch %(91.667,165)
        refiney_tiling(i,1)=refiney_1(i);
    else  %(165,200)
        refiney_tiling(i,1)=refiney1(i);

    end  
    if mod(refiney_breakpoint(i),pitch)>30/360*pitch & mod(refiney_breakpoint(i),pitch)<150/360*pitch %(18.33,91.667)
        refinex_tiling(i,1)=refinex0(i);
    elseif mod(refiney_breakpoint(i),pitch)>=150/360*pitch & mod(refiney_breakpoint(i),pitch)<270/360*pitch %(91.667,165)
        refinex_tiling(i,1)=refinex_1(i);
    else  %(165,220)
        refinex_tiling(i,1)=refinex1(i);
    end
end

uncertainty_xy=sqrt((fit_sum_Del(:,14).^2+fit_sum_Del(:,15).^2)./2)*pixelsize;
sigma=sqrt((fit_sum_Del(:,11).^2+fit_sum_Del(:,12).^2)./2)*pixelsize;
table_col={'frame','x_nm','y_nm','intensity_photon','sigma_nm','offset_photon','uncertainty_xy_nm'};
result_output=table(fit_sum_Del(:,9),refinex_tiling,refiney_tiling,fit_sum_Del(:,10),sigma,...
    fit_sum_Del(:,13),uncertainty_xy,'VariableNames',table_col) ; % The output
                                                                 % for this 
                                                                 % operation includes
                                                                 % (x,y,z,photons,frame)
                                                                 % information.    
% filename=char(filename_wformat);
% writetable(result_output,[filename(1:end-5) '_refine_sameksameb.csv']);


figure;histogram1 = histogram(refiney_tiling,'BinWidth',1);hold on
histogram(roughy,'BinWidth',1)

figure;histogram2 = histogram(refinex_tiling,'BinWidth',1);hold on
histogram(roughx,'BinWidth',1)




































