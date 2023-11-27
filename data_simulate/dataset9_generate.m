% 读取picasso-simulate生成的.raw文件
% by cxj 
% 2023.4.17
clc
clear
close all
tic
% 读取图片
load('PSF.mat');
% load('C:\Users\15600\Desktop\lattice_9ill.mat');
imgfile="C:\Users\15600\Desktop\grid10nm\";
fileID=fopen("C:\Users\15600\Desktop\grid10nm\g10nm_2000px.raw",'r');

%读取simulate生成的文件，转为.mat
g10nm_2000px_int=fread(fileID,Inf,'uint16');
fclose(fileID);
width=2000;
height=2000;
frame=numel(g10nm_2000px_int)/(width*height);
g10nm_2000px_int=reshape(g10nm_2000px_int,[width,height,frame]);
g10nm_2000px_int=permute(g10nm_2000px_int,[2,1,3]);

% 将每个800*800矩阵中的非零点数量置为1
for i = 1:frame
    % 获取当前矩阵
    currentMatrix = g10nm_2000px_int(:, :, i);
    % 找到非零点的索引位置
    [row, col] = find(currentMatrix ~= 0);
    if isempty(row)
        g10nm_2000px(:, :, i) = currentMatrix;
        continue
    end
    % 随机选择一个非零点索引位置
    randomIndex = randi(numel(row));
    % 将非零点置零
    currentMatrix(currentMatrix ~= 0) = 0;
    % 将随机选择的非零点位置设置为原来的值
    % currentMatrix(row(randomIndex), col(randomIndex)) = g10nm_2000px_int(row(randomIndex), col(randomIndex), i);
    currentMatrix(row(randomIndex), col(randomIndex)) = 1;
    % 更新当前矩阵
    g10nm_2000px(:, :, i) = currentMatrix;
end

toc
pixelsize=1;ROI=2000;pitch=330;int=5000;%都是 nm为单位，ROI/pitch~6个周期
lattice_9ill=lattice_ill(pixelsize,ROI,pitch)*int;
PSF=image_upsample(PSF,61*10);%pixelsize=1
PSF=PSF/sum(PSF(:));
gain=1000;
bk=50;
% 生成九套数据，命名为A_lattice_downsample_noise_1,2...9
% 1=pp,2=p0,3=pn,4=0p,5=00,6=0n,7=np,8=n0,9=nn
for i=1:size(lattice_9ill,3)
    tic
    for j=1:frame
        g10nm_2000px_lattice(:,:,j)=conv2(g10nm_2000px(1:1901,1:1901,j).*lattice_9ill(1:1901,1:1901,i)/9,PSF,'same');
        g10nm_2000px_lattice_downsample(:,:,j)=image_upsample(g10nm_2000px_lattice(:,:,j),20);
        %加高斯噪声
        g10nm_2000px_lattice_downsample(:,:,j)=g10nm_2000px_lattice_downsample(:,:,j);
        mean_val=max(max(g10nm_2000px_lattice_downsample(:,:,j)))/5;
        variance=mean_val/2;
        noise=normrnd(mean_val,variance,[20,20]);
        % noise=zeros(14,14);
        g10nm_2000px_lattice_downsample_noise(:,:,j)=noise+g10nm_2000px_lattice_downsample(:,:,j);
        g10nm_2000px_lattice_downsample_noise_gain(:,:,j)=g10nm_2000px_lattice_downsample_noise(:,:,j)*gain;
        g10nm_2000px_lattice_downsample_noise_gain(:,:,j)= g10nm_2000px_lattice_downsample_noise_gain(:,:,j)+bk;
        g10nm_2000px_lattice_downsample_noise_gain(:,:,j)=g10nm_2000px_lattice_downsample_noise_gain(:,:,j);
    end
    toc
    name=strcat("g10nm_2000px_bk500_",num2str(i));
    tic
    save(strcat(imgfile,name),'g10nm_2000px_lattice_downsample_noise_gain');
    toc
    % tiffwrite(g10nm_2000px_lattice_downsample_noise_gain,sprintf('%s%s.tif',imgfile,name));

    if i==1
        g10nm_2000px_lattice_downsample_noise_total=g10nm_2000px_lattice_downsample_noise_gain;
    else
        g10nm_2000px_lattice_downsample_noise_total=g10nm_2000px_lattice_downsample_noise_gain+g10nm_2000px_lattice_downsample_noise_total;
    end
end
save(strcat(imgfile,'g10nm_2000px_lattice_downsample_noise_total'),'g10nm_2000px_lattice_downsample_noise_total');

tiffwrite(g10nm_2000px_lattice_downsample_noise_total,sprintf('%s%s.tif',imgfile,'g10nm_2000px_bk500_total'));

