function higher_fidelity_simulations(choose_fig)

folder_fn = dir('*/fn_ESN/'); folder_fn = folder_fn(1).folder;
folder = [folder_fn, '/ESN_bias_data'];
%% ======================== SELECT WORKING MODELS ======================= %
model_true      =   'Wave_high-G_eqn';
model_forecast	=   'Galerkin_low-sqrt';
% --------------------------- Compare models ---------------------------- %
% Higher fidelity model (truth)
name  =   [folder,'/Truth_', model_true];
try load(name,'Sim_t'); fprintf(['Loading true data: ',model_true,'\n'])
catch; fprintf(['Creating true data: ',model_true,'\n'])
    Sim_t   =   fn_create_observations(model_true); save(name, 'Sim_t')
end
% Forecast model comparison with same parameters
name  =   [folder,'/Truth_', model_forecast];
try load(name,'Sim_f'); fprintf(['Loading forecast data: ',model_forecast,'\n'])
catch; fprintf(['Creating forecast solution: ',model_forecast,'\n'])
    Sim_f   =   fn_create_observations(model_forecast); save(name, 'Sim_f')
end
% plot_xf(Sim_t); plot_xf(Sim_f) % Plot models
% ------------------------ Compute and plot bias ------------------------ %
% specify the filename as third argument to save the data
name =	[model_true,'_+_', model_forecast];
% [B, t_B] = fn_plot_bias(Sim_t, Sim_f, [folder,'/BIAS_',name]); 

%% =========================== LOAD ESN DATA ============================ %
filename_bias   =   [folder,'\ESN_',name];
try ESN_P = load(filename_bias); fprintf('Loading ESN data \n')
catch; fprintf('Creating ESN data...\n')
    pyrunfile('ESN_RVC_validation.py', filename=name, folder=folder_fn); 
    ESN_P   =   load(filename_bias); fprintf('done!\n'); 
end

% % ----------------- Check ESN comparing to actual bias ------------------ %
% fn_check_esn(Sim_t, Sim_f, B, t_B, ESN_P)
% pause(.5);
%% ============================ ASSIMILATION ============================ %
rng(2)  % For reproducibility
% ------------------------- Initialise filter --------------------------- %
Filter          =   fn_init_filter(Sim_f, model_forecast, model_true, choose_fig);


% ------------------- Define the time of observations ------------------- %
[obs_idx, Filter] = fn_define_observation_indices(Filter, ESN_P);

% ----------------------- Apply Data Assimilation ----------------------- %
Truth.t_mic     =   Sim_t.t_mic(obs_idx);
Truth.p_mic     =   Sim_t.p_mic(:,obs_idx)';

print_assimilation_type(Filter,Sim_f, ESN_P)

if contains(model_forecast, 'Galerkin')
    [t_KF,Aa_KF,Obs,t_U,U,to,Uo] = EnSRKF_Galerkin(Truth, Filter, ESN_P);
else 
    error('EnSRKF_Wave not yet developed. Select Galerkin forecast model')
%     [t_KF, Aa_KF, Obs, t_bias, bias] = EnSRKF_Wave(Truth, Filter, ESN_P);
end
%% ============================ PLOT RESULTS ============================ %
fn_plot_results(Filter, Sim_t, Truth, t_KF, Aa_KF, Obs, t_U, U, to, Uo)

end


