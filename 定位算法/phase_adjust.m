function fit_phase_adjust=phase_adjust(photon_phase)
    if any(photon_phase<0)
        if norm(max(photon_phase)-min(photon_phase))>0.5712
            neg_id=find(photon_phase<0); 
            fit_phase_adjust=photon_phase;
            fit_phase_adjust(neg_id)=fit_phase_adjust(neg_id)+2*pi;%令负值取正
        
        elseif sum(photon_phase<-0.2)>1
            fit_phase_adjust=photon_phase+ones(1,3)*2*pi;
        else
            fit_phase_adjust=photon_phase;
        end
    else
        fit_phase_adjust=photon_phase;
    end
end