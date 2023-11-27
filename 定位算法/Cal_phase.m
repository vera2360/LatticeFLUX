function fit_phase=Cal_phase(fit_sum)
    imglen=size(fit_sum,1);
    for n=1:imglen

        [py_p,py_0,py_n,mdy_p,mdy_0,mdy_n]=cal_phasey(fit_sum(n,5:13));%p代表+2/3*pi，0代表0，-代表-2/3*pi
        [px_p,px_0,px_n,mdx_p,mdx_0,mdx_n]=cal_phasex(fit_sum(n,5:13));
        fit_phase(n,:)=[py_p,py_0,py_n,px_p,px_0,px_n,mdy_p,mdy_0,mdy_n,mdx_p,mdx_0,mdx_n];
        
    end
end