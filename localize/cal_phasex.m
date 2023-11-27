function [px1,px2,px3,mdx1,mdx2,mdx3]=cal_phasex(phasex_par)
    %注意147一组，258一组，369一组
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
%% 计算px1，px1是不受系数1+mdy*sin(py)影响，此时y相位为+2/3*pi
    offsetx1 = (phot4*cp - (phot1+phot7)/2)/(cp-1);
    phot1 = phot1 - offsetx1;
    phot4 = phot4 - offsetx1;
    phot7 = phot7 - offsetx1;
    
    ampx1 = sqrt(phot4^2 + (phot7-phot4*cp)^2/(sp^2));
    
    sx1 = phot4/ampx1;
    cx1 = (-phot7 + phot4*cp)/sp/ampx1;
    
    tpx1 = acos(cx1);
    if sx1 <0
        px1 = -real(tpx1);
    else
        px1 = real(tpx1);
    end
%% 计算px2,px2同样不受系数1+mdy*sin(py)的影响,此时y相位为0
    offsetx2 = (phot5*cp - (phot2+phot8)/2)/(cp-1);
    phot2 = phot2 - offsetx2;
    phot5 = phot5 - offsetx2;
    phot8 = phot8 - offsetx2;
    
    ampx2 = sqrt(phot5^2 + (phot8-phot5*cp)^2/(sp^2));
    
    sx2 = phot5/ampx2;
    cx2 = (-phot8 + phot5*cp)/sp/ampx2;
    
    tpx2 = acos(cx2);%cy在[-1,1]之间，计算结果为[0,pi]
    if sx2 <0
        px2 = -real(tpx2); %sin<0,表示在[-pi,0]之间，但通过轴对称特性，acos的计算值可以直接扩展到[-pi,pi]
    else
        px2 = real(tpx2);
    end
%% 计算px3,px3同样不受系数1+mdy*sin(py)的影响,此时y相位为-2/3*pi
    offsetx3 = (phot6*cp - (phot3+phot9)/2)/(cp-1);
    phot3 = phot3 - offsetx3;
    phot6 = phot6 - offsetx3;
    phot9 = phot9 - offsetx3;
    
    ampx3 = sqrt(phot6^2 + (phot9-phot6*cp)^2/(sp^2));
    
    sx3 = phot6/ampx3;
    cx3 = (-phot9 + phot6*cp)/sp/ampx3;
    
    tpx3 = acos(cx3);%cy在[-1,1]之间，计算结果为[0,pi]
    if sx3 <0
        px3 = -real(tpx3); %sin<0,表示在[-pi,0]之间，但通过轴对称特性，acos的计算值可以直接扩展到[-pi,pi]
    else
        px3 = real(tpx3);
    end
    %% mdx和mdy通过amp/offset可以求的，phasex和phasey的求解也不受系数的影响
    mdx1=ampx1/offsetx1;
    mdx2=ampx2/offsetx2;
    mdx3=ampx3/offsetx3;
end