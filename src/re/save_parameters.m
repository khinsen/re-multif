#!/usr/bin/octave -qf
% usage: octave save_parameters.m

% experimental context parameters
experiment_class = "period"; % "freqdiv" | "period"
period_experiment_range = [2e4, 0.5e4, 16e4]; % min : step : max

% initial concentrations of repressor proteins
R1 = 50e-9; % (M)
R2 = 50e-9; % (M)
R3 = 0e-9; % (M)
R4 = 0e-9; % (M)

% input signal configuration
input_type = "sine"; % "square" | "sine"
input_dc_level = 6e-9; % (M)
input_period = 0.9e5; % (s) can be overrided with input_period_range
input_amplitude = 50e-9; % (M)
input_duty_cycle = 0.5; % (normalized %)

% simulation time controls
simulation_total_time = 4.0e5; %(s)
simulation_step = 60; %(s)

% underlying model type
model = "Full"; % "QSSA" | "Full"

% plot configuration
plot_output_filename = "freqdiv-period.pdf"; % "filename.extension"
plot_include_input_signal = false; % (bool)

% network parameters (follows the nomenclature in the original paper)
k_tl = 6e-4; % Translation rate (1/s)
delta_m = 2.5e-3; % mRNA degradation rate (1/s)
delta_x = 4e-4; % Protein degradation rate (1/s)
N = 1.3; % Hill coefficient (scalar)
k_a = 2e-8; % k_a for activator (M)
k_r = 6e-10; % k_a for repressor (M)
beta = 4e-10; % Maximum transcription rate (P1 and P2)
P_tc = 4e-10; % Unrepressed transcription rate (P3–P6)

% rename output'ed "parameters.mat" to something more descriptive of each sim.
save "parameters.mat" ...
     experiment_class period_experiment_range ...
     R1 R2 R3 R4 ...
     input_type input_dc_level input_period input_amplitude input_duty_cycle ...
     simulation_total_time simulation_step ...
     model ...
     plot_output_filename plot_include_input_signal ...
     k_tl delta_m delta_x N k_a k_r beta P_tc;
