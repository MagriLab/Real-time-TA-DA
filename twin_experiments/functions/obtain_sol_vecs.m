function V = obtain_sol_vecs(SOLS,P)
% Simple function that computes the pressure and velocity at the flame
% location from the solution vectors from Function_main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



N_m         =   P.N_m;
sin_jpixf   =   P.sin_jpixf;
cos_jpixf   =   P.cos_jpixf;

%%%%%%%%%%%%%%%%%%%%%%%%%%% TRUE STATE SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
if  isfield(SOLS, 'Aa_true') && ~isempty(SOLS.Aa_true) 
    Aa_true 	=   SOLS.Aa_true;
    eta_true    =   Aa_true(1:N_m,:);
    etad_true   =   Aa_true(N_m+1:2*N_m,:);

    v_xf_true	=   cos_jpixf * eta_true;    
    p_xf_true	=   - sin_jpixf * etad_true;

    V.eta_true	=   eta_true;
    V.etad_true	=   etad_true;
    V.v_xf_true	=   v_xf_true;
    V.p_xf_true	=   p_xf_true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NO DA SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  isfield(SOLS, 'Aa_comp') && ~isempty(SOLS.Aa_comp) 
    Aa_comp     =   SOLS.Aa_comp;
    eta_comp    =   Aa_comp(1:N_m,:);
    etad_comp   =   Aa_comp(N_m+1:2*N_m,:);

    v_xf_comp	=   cos_jpixf * eta_comp;    
    p_xf_comp	=   - sin_jpixf * etad_comp;

    V.eta_comp	=   eta_comp;
    V.etad_comp	=   etad_comp;
    V.v_xf_comp	=   v_xf_comp;
    V.p_xf_comp	=   p_xf_comp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% FILTERED SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
if  isfield(SOLS, 'Aa')
    Aa      =   SOLS.Aa;
    eta     =   mean(Aa(:,1:N_m,:),3);       % Mean in the 3-direction i.e. n_ens
    etad    =   mean(Aa(:,N_m+1:2*N_m,:),3); % Mean in the 3-direction i.e. n_ens

    v_xf  	=   cos_jpixf * eta';    
    p_xf  	=   - sin_jpixf * etad';

    V.eta  	=   eta;
    V.etad	=   etad;
    V.v_xf 	=   v_xf;
    V.p_xf 	=   p_xf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER ESTIAMATION %%%%%%%%%%%%%%%%%%%%%%%%%%
if P.param_est == 2
    V.tau_PE	=   mean(Aa(:,end,:),3);
    V.beta_PE	=   mean(Aa(:,end-1,:),3); 
elseif P.param_est == 1
    V.beta_PE	=   mean(Aa(:,end,:),3); 
end

end

