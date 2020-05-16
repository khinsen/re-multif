#!/usr/bin/octave -qf
% usage: octave run_experiments.m path/to/parameters.m

clear;
format long;


% load model parameters
parameter_file = argv(){1};
load(parameter_file);

% calculate QSSA equation coefficients
alpha = k_tl * beta / delta_m;
gamma = k_tl * P_tc / delta_m;

% short-hand for Hill functions
Ha = @(X) hill_activation(X, N, k_a);
Hr = @(X) hill_repression(X, N, k_r);

% choose input signal
if strcmp(input_type, "square")
  I = @(t) square_wave(t, input_period, input_amplitude, ...
                       input_dc_level, input_duty_cycle);
else
  error("Undefined input type: %s.", input_type)
endif

% choose the network's internal model
if strcmp(model, "QSSA")
  model = @(R, t) qssa(R, t, I, Ha, Hr, alpha, gamma, delta_x);
else
  error("Undefined model type: %s.", model)
endif

% set initial state
R = [R1, R2, R3, R4];

% run the simulation
step_number = floor(simulation_total_time / simulation_step) + 1;
Rs = euler_simulate(model, R, step_number, simulation_step);


figure;

% scale data for easier visualization
timescale = 1e5; % 10^5 seconds
molscale = 1e-9; % nM

% time, repressor concentration and input concentration series
t = (0 : simulation_step : simulation_total_time)' / timescale;
R1s = Rs(:,1) / molscale;
R2s = Rs(:,2) / molscale;
R3s = Rs(:,3) / molscale;
R4s = Rs(:,4) / molscale;
Is = zeros(size(t));
for i = 1 : length(t)
  Is(i) = I(i * simulation_step) / molscale;
endfor

% plot input behaviour (if needed)
if plot_include_input_signal
  subplot(2, 1, 2);
  plot(t, Is, 'b;I;');
  axis([-Inf, +Inf, 0, input_amplitude * 1.23 / molscale]);
  pbaspect([1 0.334 1]);
  xlabel("Time (10^5 seconds)");
  ylabel("Concentration (nM)");

  subplot(2, 1, 1); % subplot() may or may not be needed for the model plot
endif

% plot model behaviour
plot(t,R1s,'--m;R1;', t,R2s,':k;R2;', t,R3s,'-r;R3;', t,R4s,'-.g;R4;');
pbaspect([1 0.334 1]);
xlabel("Time (10^5 seconds)");
ylabel("Concentration (nM)");

% save plotted figure (gets put in the same folder the script runs from)
print(plot_output_filename);
