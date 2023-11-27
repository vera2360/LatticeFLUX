function fit_locs=fit_9subset_MLE(imgbuf,cx,cy,sigmax,sigmay,r,frame)
    [xx,yy]=meshgrid(round(cx-r):round(cx+r),round(cy-r):round(cy+r));
    x=double([xx,yy]);
    ff=double(imgbuf);
    bklist=min(min(ff,[],1),[],2);
    bklist=bklist(:);
    peaklist=max(max(ff,[],1),[],2);
    peaklist=peaklist(:);
    xstart=double([cx,cy,sigmax,sigmay,peaklist',bklist']);
    options=optimset('Display','off','MaxFunEvals',1e7,'MaxIter',100,'TolFun',0.0001,'Algorithm','trust-region-reflective');
    [fit_locs,~,exitflag]=MLEcurvefit(@Gaussian6,xstart,x,ff,options);


end


function f=Gaussian6(xstart,x)
    xd=x(:,1:9);
    yd=x(:,10:18);
    f1=xstart(5)*(exp(-0.5*(xd-xstart(1)).^2./(xstart(3)^2)-0.5*(yd-xstart(2)).^2./(xstart(4)^2)))+xstart(14);
    f2=xstart(6)*(exp(-0.5*(xd-xstart(1)).^2./(xstart(3)^2)-0.5*(yd-xstart(2)).^2./(xstart(4)^2)))+xstart(15);
    f3=xstart(7)*(exp(-0.5*(xd-xstart(1)).^2./(xstart(3)^2)-0.5*(yd-xstart(2)).^2./(xstart(4)^2)))+xstart(16);
    f4=xstart(8)*(exp(-0.5*(xd-xstart(1)).^2./(xstart(3)^2)-0.5*(yd-xstart(2)).^2./(xstart(4)^2)))+xstart(17);
    f5=xstart(9)*(exp(-0.5*(xd-xstart(1)).^2./(xstart(3)^2)-0.5*(yd-xstart(2)).^2./(xstart(4)^2)))+xstart(18);
    f6=xstart(10)*(exp(-0.5*(xd-xstart(1)).^2./(xstart(3)^2)-0.5*(yd-xstart(2)).^2./(xstart(4)^2)))+xstart(19);
    f7=xstart(11)*(exp(-0.5*(xd-xstart(1)).^2./(xstart(3)^2)-0.5*(yd-xstart(2)).^2./(xstart(4)^2)))+xstart(20);
    f8=xstart(12)*(exp(-0.5*(xd-xstart(1)).^2./(xstart(3)^2)-0.5*(yd-xstart(2)).^2./(xstart(4)^2)))+xstart(21);
    f9=xstart(13)*(exp(-0.5*(xd-xstart(1)).^2./(xstart(3)^2)-0.5*(yd-xstart(2)).^2./(xstart(4)^2)))+xstart(22);
    f = cat(3, f1,f2,f3,f4,f5,f6,f7,f8,f9);
end