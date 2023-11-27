function lattice_9ill=lattice_ill(pixelsize,ROI,pitch)
%%%%% ���������ⳡ������330 nm���ⳡ�ܳߴ�2560 nm�����سߴ�1 nm

    I0=1;%�ⳡƽ����ǿ
    m=1;%�������
    x=0:pixelsize:(ROI-1);y=0:pixelsize:(ROI-1);%�����������ص�ߴ�
    [X,Y]=meshgrid(x,y);
    
    %% ���ɾ���������ⳡ
    px1 = 0; py1 = 0;%��ʼ��λ
    lattice_00 = I0*(1+m*sin(2*pi/pitch*X+px1)).*(1+m*sin(2*pi/pitch*Y+py1));

    px2 = pi*2/3; py1 = 0;%��ʼ��λ
    lattice_p0 = I0*(1+m*sin(2*pi/pitch*X+px2)).*(1+m*sin(2*pi/pitch*Y+py1));

    px3 = pi*4/3; py1 = 0;%��ʼ��λ
    lattice_n0 = I0*(1+m*sin(2*pi/pitch*X+px3)).*(1+m*sin(2*pi/pitch*Y+py1));

    px1 = 0; py2 = pi*2/3;%��ʼ��λ
    lattice_0p = I0*(1+m*sin(2*pi/pitch*X+px1)).*(1+m*sin(2*pi/pitch*Y+py2));

    px2 = pi*2/3; py2 = pi*2/3;%��ʼ��λ
    lattice_pp = I0*(1+m*sin(2*pi/pitch*X+px2)).*(1+m*sin(2*pi/pitch*Y+py2));

    px3 = pi*4/3; py2 = pi*2/3;%��ʼ��λ
    lattice_np = I0*(1+m*sin(2*pi/pitch*X+px3)).*(1+m*sin(2*pi/pitch*Y+py2));

    px1 = 0;py3 = pi*4/3;%��ʼ��λ
    lattice_0n = I0*(1+m*sin(2*pi/pitch*X+px1)).*(1+m*sin(2*pi/pitch*Y+py3));

    px2 = pi*2/3; py3 = pi*4/3;%��ʼ��λ
    lattice_pn = I0*(1+m*sin(2*pi/pitch*X+px2)).*(1+m*sin(2*pi/pitch*Y+py3));

    px3 = pi*4/3;py3 = pi*4/3;%��ʼ��λ
    lattice_nn = I0*(1+m*sin(2*pi/pitch*X+px3)).*(1+m*sin(2*pi/pitch*Y+py3));
    %ע�����ӵ�˳������ɵ������˳��ͬ
    lattice_9ill=cat(3,lattice_pp,lattice_p0,lattice_pn,lattice_0p,lattice_00,lattice_0n,lattice_np,lattice_n0,lattice_nn);
   
%     subplot(3,3,1);pcolor(lattice_00(:,:));shading flat;axis equal off;colormap(jet);
%     subplot(3,3,2);pcolor(lattice_p0(:,:));shading flat;axis equal off;colormap(jet);
%     subplot(3,3,3);pcolor(lattice_n0(:,:));shading flat;axis equal off;colormap(jet);
%     subplot(3,3,4);pcolor(lattice_0p(:,:));shading flat;axis equal off;colormap(jet);
%     subplot(3,3,5);pcolor(lattice_pp(:,:));shading flat;axis equal off;colormap(jet);
%     subplot(3,3,6);pcolor(lattice_np(:,:));shading flat;axis equal off;colormap(jet);
%     subplot(3,3,7);pcolor(lattice_0n(:,:));shading flat;axis equal off;colormap(jet);
%     subplot(3,3,8);pcolor(lattice_pn(:,:));shading flat;axis equal off;colormap(jet);
%     subplot(3,3,9);pcolor(lattice_nn(:,:));shading flat;axis equal off;colormap(jet);
    %% �洢����
    save('lattice_9ill.mat','lattice_9ill');
end

