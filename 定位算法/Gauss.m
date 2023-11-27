% function result=Gauss(pixel)
%pixel决定了光斑大小，origin=50
clc;clear;
    syms theta phi;
    n=1.515;
    NA=1.49;
    lambda=680;%um,emission wavelength
    k=2.*pi/lambda;
    z2=0;%焦点处

    pixelsize=10;%nm
    x=-300:pixelsize:300;
    x(x==0)=eps;
    y=-300:pixelsize:300;

    lenx=length(x);leny=length(y);

    Ex=zeros(lenx,leny); Ey=zeros(lenx,leny);
    % matlabpool open;
    tic
    for ii=1:lenx
        parfor jj=1:leny
            r2=sqrt(x(ii).^2+y(jj).^2);
            phi2=atan(y(jj)./x(ii));
            %% 高斯光斑
             FEx=@(theta,phi)sin(theta).*sqrt(cos(theta)).*((2^(1/2)*((cos(theta) - 1)*cos(phi)^2 + 1))/2 + (2^(1/2)*cos(phi)*sin(phi)*(cos(theta) - 1)*i)/2).*exp(i.*k.*n.*(z2.*cos(theta)+r2.*sin(theta).*cos(phi-phi2)));
             FEy=@(theta,phi)sin(theta).*sqrt(cos(theta)).*((2^(1/2)*((cos(theta) - 1)*sin(phi)^2 + 1)*i)/2 + (2^(1/2)*cos(phi)*sin(phi)*(cos(theta) - 1))/2).*exp(i.*k.*n.*(z2.*cos(theta)+r2.*sin(theta).*cos(phi-phi2)));
             FEz=@(theta,phi)sin(theta).*sqrt(cos(theta)).*((2^(1/2)*cos(phi)*sin(theta))/2 + (2^(1/2)*sin(phi)*sin(theta)*i)/2).*exp(i.*k.*n.*(z2.*cos(theta)+r2.*sin(theta).*cos(phi-phi2)));         
             Ex(ii,jj)=dblquad(FEx,0,asin(NA./n),0,2.*pi);
             Ey(ii,jj)=dblquad(FEy,0,asin(NA./n),0,2.*pi);
             Ez(ii,jj)=dblquad(FEz,0,asin(NA./n),0,2.*pi);
        end
    end
    toc
    % matlabpool close;
    E=(abs(i.*Ex)).^2+(abs(i.*Ey)).^2+(abs(i.*Ez)).^2;
    result=E./max(max(E));
    x=-300:pixelsize:300;y=-300:pixelsize:300;
    pcolor(y,x,result); %绘制伪彩图
    colormap(gray);
    shading flat; %去掉图上的网格线
    axis equal ;%等比坐标轴
%end
