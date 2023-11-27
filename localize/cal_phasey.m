function [py1,py2,py3,mdy1,mdy2,mdy3]=cal_phasey(phasex_par)
%注意求y相位，会求出3个来，一行求出一个
    phot1=phasex_par(1);
    phot2=phasex_par(2);
    phot3=phasex_par(3);
    phot4=phasex_par(4);
    phot5=phasex_par(5);
    phot6=phasex_par(6);
    phot7=phasex_par(7);
    phot8=phasex_par(8);
    phot9=phasex_par(9);

    cp=cos(120/180*pi);
    sp=sin(120/180*pi);% 已知phase shift=120

%% 计算py1，py1是不受系数1+mdx*sin(px)影响，此时x相位为+2/3*pi
    offsety1 = (phot2*cp - (phot1+phot3)/2)/(cp-1);
    phot1 = phot1 - offsety1;
    phot2 = phot2 - offsety1;
    phot3 = phot3 - offsety1;
    
    ampy1 = sqrt(phot2^2 + (phot3-phot2*cp)^2/(sp^2));
    
    sy1 = phot2/ampy1;
    cy1 = (-phot3 + phot2*cp)/sp/ampy1;
    
    tpy1 = acos(cy1);
    if sy1 <0
        py1 = -real(tpy1);
    else
        py1 = real(tpy1);
    end
%% 计算py2,py2同样不受系数1+mdx*sin(px)的影响,此时x相位为0
    offsety2 = (phot5*cp - (phot4+phot6)/2)/(cp-1);
    phot4 = phot4 - offsety2;
    phot5 = phot5 - offsety2;
    phot6 = phot6 - offsety2;
    
    ampy2 = sqrt(phot5^2 + (phot6-phot5*cp)^2/(sp^2));
    
    sy2 = phot5/ampy2;
    cy2 = (- phot6 + phot5*cp)/sp/ampy2;
    
    tpy2 = acos(cy2);%cy在[-1,1]之间，计算结果为[0,pi]
    if sy2 <0
        py2 = -real(tpy2); %sin<0,表示在[-pi,0]之间，但通过轴对称特性，acos的计算值可以直接扩展到[-pi,pi]
    else
        py2 = real(tpy2);
    end
%% 计算py3,py3同样不受系数1+mdx*sin(px)的影响,此时x相位为-2/3*pi
    offsety3 = (phot8*cp - (phot7+phot9)/2)/(cp-1);
    phot7 = phot7 - offsety3;
    phot8 = phot8 - offsety3;
    phot9 = phot9 - offsety3;
    
    ampy3 = sqrt(phot8^2 + (phot9-phot8*cp)^2/(sp^2));
    
    sy3 = phot8/ampy3;
    cy3 = (-phot9 + phot8*cp)/sp/ampy3;
    
    tpy3 = acos(cy3);%cy在[-1,1]之间，计算结果为[0,pi]
    if sy3 <0
        py3 = -real(tpy3); %sin<0,表示在[-pi,0]之间，但通过轴对称特性，acos的计算值可以直接扩展到[-pi,pi]
    else
        py3 = real(tpy3);
    end
    %% mdx和mdy通过amp/offset可以求的，phasex和phasey的求解也不受系数的影响
    mdy1=ampy1/offsety1;
    mdy2=ampy2/offsety2;
    mdy3=ampy3/offsety3;

    %% 先暂时不求调制深度、幅值、偏置
%     ay=ampy/(1+mdx*sin(px));
%     oy=offsety/(1+mdx*sin(px));
% 
%     ax=ampy/(1+mdy*sin(py));
%     ox=ampy/(1+mdy*sin(py));
end