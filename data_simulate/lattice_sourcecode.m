clear
e_1 = 1;e_2 =1;e_3 =1;e_4 =1;
lamda=488;
NA=1.49;
n=1.517;
k=2*pi/lamda; 
theta=asin(NA/n);


deltax=0;deltay=0;deltaz=0;

pixel=200; 
dx=200./pixel;% um  1nm
xx=0:dx:1400;
yy=0:dx:1400;
zz=0:dx:1400;
lenx=length(xx); 
leny=length(yy);
lenz=length(zz);
CoreNum=10;
% if isempty(gcp('nocreat'))
%     parpool(CoreNum);
% end
for kk=5:lenz
    for jj=1:leny
        for ii=1:lenx
            I(ii,jj,kk)=(e_1.^2+e_2.^2+e_3.^2+e_4.^2)+2*(e_1*e_2*cos(0)*cos(-2*k*sin(theta)*(xx(ii)+deltax))...
                +e_1*e_3*cos(0)*cos(-k*sin(theta)*(xx(ii)+deltax)+k*sin(theta)*(yy(jj)+deltay))+e_1*e_4*cos(0)*cos(-k*sin(theta)*(xx(ii)+deltax)-k*sin(theta)*(yy(jj)+deltay))...
                +e_2*e_3*cos(0)*cos(k*sin(theta)*(xx(ii)+deltax)+k*sin(theta)*(yy(jj)+deltay))...
                +e_2*e_4*cos(0)*cos(k*sin(theta)*(xx(ii)+deltax)-k*sin(theta)*(yy(jj)+deltay))+e_3*e_4*cos(0)*cos(-2*k*sin(theta)*(yy(jj)+deltay)));

        end
    end
    T=I(:,:,kk);subplot(121);imagesc(T);axis equal;shading flat
end



            
            
            
            