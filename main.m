clc; addpath(genpath(pwd)); plot_settings
% ----------------------------------------------------------------------- %

choose_fig  =   24;      % Select any of the figures in the paper

% ----------------------------------------------------------------------- %
if choose_fig <= 22
    twin_experiments(choose_fig)
else
    higher_fidelity_simulations(choose_fig)
end
% ----------------------------------------------------------------------- %