#!/usr/bin/octave -qf
% usage: octave save_parameters.m
% then rename the generated "parameters.mat" file

% experimental context parameters
experiment_class = "switch"; % "single" | "period" | "oscillator | "switch"
model = "Full"; % "QSSA" | "Full"
period_experiment_range = [2e4, 0.5e4, 16e4]; % min : step : max
oscillator_input_range = [0.4e-9, 0.1e-9, 7e-9]; % min : step : max
switch_trigger_delta = +399941e-10; % k_r
switch_trigger_time = [1.50e5 1.55e5]; % (s) start, stop
switch_R = [1 1 0 0]; % affects R1? R2? R3? R4?

% initial concentrations of repressor proteins, overrided with bifurcation_init
R1 = 50e-9; % (M)
R2 = 50e-9; % (M)
R3 = 0e-9; % (M)
R4 = 0e-9; % (M)

% input signal configuration
input_type = "constant"; % "sine" | "square" | "constant"
input_dc_level = 0.1e-9; % (M), may be overrided with oscillator_input_range
input_period = 0.9e5; % (s), may be overrided with input_period_range
input_amplitude = 50e-9; % (M)
input_duty_cycle = 0.5; % (normalized %)

% simulation time controls
simulation_total_time = 3.5e5; %(s)
simulation_step = 60; %(s)

% plot configuration
plot_output_filename = "switch-A.pdf"; % `filename.extension`
plot_include_input_signal = false; % (bool)

% save experiment configuration in plain text
save "parameters.mat" ...
     experiment_class model ...
     period_experiment_range ...
     oscillator_input_range ...
     switch_trigger_delta switch_trigger_time switch_R ...
     R1 R2 R3 R4 ...
     input_type input_dc_level input_period input_amplitude input_duty_cycle ...
     simulation_total_time simulation_step ...
     plot_output_filename plot_include_input_signal;
