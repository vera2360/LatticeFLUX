function fit_phase_period=Cal_period(gauss_position,pitch,pixelsize,phase)

    period=gauss_position*pixelsize/pitch;
    photon_phase=phase+2*pi*floor(period)*ones(1,3);
    gauss_phase=gauss_position.*pixelsize./pitch.*2*pi;%一开始的时候没有将高斯的定位坐标换算为相位坐标
    fit_phase_period=photon_phase;
    for i=1:3
        if gauss_phase-photon_phase(i)>3.5
            fit_phase_period(i)=photon_phase(i)+2*pi;
        elseif gauss_phase-photon_phase(i)<-3.5
            fit_phase_period(i)=photon_phase(i)-2*pi;
        else
            fit_phase_period(i)=photon_phase(i);
        end
    end
end