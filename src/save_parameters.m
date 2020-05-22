#!/usr/bin/octave -qf
% usage: octave save_parameters.m
% then rename the generated "parameters.mat" file

% experimental context parameters
experiment_class = "oscillator"; % "single" | "period" | "oscillator | "switch"
model = "Full"; % "QSSA" | "Full"
period_experiment_range = [2e4, 0.5e4, 16e4]; % min : step : max
oscillator_input_range = [0.4e-9, 0.1e-9, 7e-9]; % min : step : max

% initial concentrations of repressor proteins, overrided with bifurcation_init
R1 = 50e-9; % (M)
R2 = 50e-9; % (M)
R3 = 0e-9; % (M)
R4 = 0e-9; % (M)

% input signal configuration
input_type = "constant"; % "sine" | "square" | "constant"
input_dc_level = 6e-9; % (M)
input_period = 0.9e5; % (s) can be overrided with input_period_range
input_amplitude = 50e-9; % (M)
input_duty_cycle = 0.5; % (normalized %)

% simulation time controls
simulation_total_time = 4.0e5; %(s)
simulation_step = 60; %(s)

% plot configuration
plot_output_filename = "freqdiv-period.pdf"; % "filename.extension"
plot_include_input_signal = false; % (bool)

% save experiment configuration in plain text
save "parameters.mat" ...
     experiment_class model period_experiment_range oscillator_input_range ...
     R1 R2 R3 R4 ...
     input_type input_dc_level input_period input_amplitude input_duty_cycle ...
     simulation_total_time simulation_step ...
     plot_output_filename plot_include_input_signal;
