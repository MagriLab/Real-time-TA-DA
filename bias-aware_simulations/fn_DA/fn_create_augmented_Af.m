function [Af_aug, Filter] = fn_create_augmented_Af(Af,pA,SE,PE,Filter)
    
    m       =   Filter.m;
    N       =   Filter.N;
    N_c     =   Filter.N_c;
    N_m     =   Filter.N_m;
    N_mic   =   Filter.N_mic;
    Ni_state	=   1:2*N_m+N_c;
    Ni_params  	=   2*N_m+N_c+1:N;
    % Make sure Af is m x N
    if size(Af,2) ~= m
        Af = Af';
    end
    % Keep state and/or parameter only in the state vector
    indices     =   {Ni_state, Ni_params};
    keep_idx    =   cell2mat(indices(find([SE,PE])));
    Af_keep     =   Af(keep_idx, :);
%     if all([SE,PE])
%         Af_keep    =	Af;
%     elseif SE
%         Af_keep    =	Af(Ni_state,:);
%     elseif PE
%         Af_keep    =	Af(Ni_params,:);
%     end
    Af_aug      = 	vertcat(Af_keep,pA);
    % Define observation matrix    
    Filter.y_0	=   horzcat(zeros(1,size(Af_keep,1)), ...
                            ones(1,Filter.N_mic));
    M           =   zeros(N_mic,length(Filter.y_0));   % 
    iq          =   1;
    for i = 1:length(Filter.y_0)
        if Filter.y_0(i) ~= 0  
            M(iq,i) =   1;      
            iq      =   iq + 1;    
        end
    end
    Filter.M    =   M;
end