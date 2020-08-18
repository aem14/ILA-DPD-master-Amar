clear; clc; close all;
%% Setup Everything
 
% Add the submodules to path
addpath(genpath('OFDM-Matlab'))
addpath(genpath('WARPLab-Matlab-Wrapper'))
addpath(genpath('Power-Amplifier-Model'))
 
rms_input = 0.50;
 
% Setup the PA simulator or TX board
PA_board = 'webRF'; % either 'WARP', 'webRF', or 'none'
 
switch PA_board
    case 'WARP'
        warp_params.nBoards = 1;         % Number of boards
        warp_params.RF_port  = 'A2B';    % Broadcast from RF A to RF B. Can also do 'B2A'
        board = WARP(warp_params);
        Fs = 40e6;    % WARP board sampling rate.
    case 'none'
        board = PowerAmplifier(7, 4);
        Fs = 40e6;    % WARP board sampling rate.
    case 'webRF'
        dbm_power = -25; % Originally -22
        board = webRF(dbm_power);
end
 
% Setup OFDM
ofdm_params.nSubcarriers = 1200;
ofdm_params.subcarrier_spacing = 15e3; % 15kHz subcarrier spacing
ofdm_params.constellation = 'QPSK';
ofdm_params.cp_length = 144; % Number of samples in cyclic prefix.
ofdm_params.nSymbols = 14;
modulator = OFDM(ofdm_params);
 
% Create TX Data
[tx_data, ~] = modulator.use;
tx_data = Signal(tx_data, modulator.sampling_rate, rms_input);
tx_data.upsample(board.sample_rate)
 
% Setup DPD
dpd_params.order = 3;
dpd_params.memory_depth = 1;
dpd_params.lag_depth = 0;  % 0 is a standard MP. >0 is GMP.
dpd_params.nIterations = 1;
dpd_params.learning_rate = 1;
dpd_params.learning_method = 'newton'; % Or 'ema' for exponential moving average.
dpd_params.use_even = false; 
dpd_params.use_conj = 0;    % Conjugate branch. Currently only set up for MP (lag = 0)
dpd_params.use_dc_term = 0; % Adds an additional term for DC
dpd = ILA_DPD(dpd_params);
 
%% Run Experiment
[~, w_out_dpd] = board.transmit(tx_data.data);
dpd.perform_learning(tx_data.data, board);
original_dpd_coeffs = dpd.coeffs;

[~, w_dpd] = board.transmit(dpd.predistort(tx_data.data));
after = w_dpd.measure_all_powers;

% Plot the before and after
w_out_dpd.name = 'Without DPD after RFWebLab';
w_out_dpd.plot_psd;

w_dpd.name = 'LS Solution DPD after RFWebLab';
w_dpd.plot_psd;


% Record the L1 power
after_values = [after(1,1)];
gradient_values = [];
hess_values = [];
total_gradient = inf;

%% Complex grid search 
step_size = 0.05; 


test_index = 1;
best_result = inf; 
best_index = 0;
coeff_array = [];
results_array = [];
for real_diff = -2*step_size:step_size:2*step_size
    for imag_diff = -2*step_size:step_size:2*step_size
        this_complex_diff = real_diff + 1i * imag_diff;
        dpd.coeffs = original_dpd_coeffs;
        dpd.coeffs(2) = dpd.coeffs(2) + this_complex_diff;
        coeff_array = [coeff_array dpd.coeffs(2)];  % Append this coeff to array
        [~, test_dpd] = board.transmit(dpd.predistort(tx_data.data));
        this_power_results = test_dpd.measure_all_powers;
        
        % Check to see if this is the new min! 
        if (this_power_results(1) < best_result)
            best_result = this_power_results(1);
            best_index = test_index;
        end
        
        % Save all the results
        results_array = [results_array this_power_results(1)];
        test_index = test_index + 1;
        save('checkpoint.mat');
    end
end

%% Plot the results
% TODO.
