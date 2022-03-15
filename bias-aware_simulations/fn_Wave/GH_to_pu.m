function [t, P_mic, U_mic] = GH_to_pu(Sim,ts,Gs,Hs)
% ----------------------------------------------------------------------- %
    Geom        = Sim.Geom;
    Mean        = Sim.Mean;
    Measurement = Sim.Measurement;
    fs          = Measurement.Fs;
    t           = ts(1):1/fs:ts(end);
    x_mic = Measurement.Mic_Pos;
    Acoustic_only = true;
    Lu  = Geom.Lu;
    Lb  = Geom.Lb;
    P_mic = zeros(length(x_mic),length(t));
   	U_mic = zeros(length(x_mic),length(t));
    for k = 1:length(x_mic)
        if x_mic(k) < -Lu || x_mic(k) > Lb
            % If location is out of the combustor, send white noise:
            P_mic(k,:) = randn(1,length(t))*10^-3;
        else
            for j = 1:length(t)
                idx = [max(j-100,1), min(j+100, length(ts))];
                [uu,pp] = fn_modehshapes(t(j),x_mic(k),...
                            ts(idx(1):idx(2)), Gs(idx(1):idx(2)), Hs(idx(1):idx(2)),...%                             ts,Gs,Hs,...
                            Mean,Geom,Acoustic_only);
                P_mic(k,j) = pp;
                U_mic(k,j) = uu;
            end
        end
    end
% ----------------------------------------------------------------------- %
end