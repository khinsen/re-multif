#!/usr/bin/octave -qf
% usage: octave run_experiments.m network.mat experiment.mat

clear;
format long;

pkg load signal; % for findpeaks


% parameterize experiment
params = argv();
load(params{1}); % network reaction constants
load(params{2}); % experiment configuration

% set initial protein concentrations
R = [R1 R2 R3 R4];

% choose input signal shape (period T left unfixed)
if strcmp(input_type, "sine")
  I = @(t, T) input_amplitude * (-cos(t * 2*pi/T) + 1) / 2 + input_dc_level;
elseif strcmp(input_type, "square")
  I = @(t, T) square_wave(t, T, ...
                          input_amplitude, input_dc_level, input_duty_cycle);
elseif strcmp(input_type, "constant")
  I = @(t, T) input_dc_level;
  input_period = Inf;
else
  error("Undefined input type: %s.", input_type);
endif

% define short-hand for Hill functions
k_a = [k_a k_a k_a k_a k_a];
k_r = [k_r k_r k_r k_r];
Ha = @(X, k) hill_activation(X, N, k);
Hr = @(X, k) hill_repression(X, N, k);

% choose the network's internal model (input I and constant k_r left unfixed)
if strcmp(model, "QSSA")
  alpha = k_tl * beta / delta_m;
  gamma = k_tl * P_tc / delta_m;
  model = @(R, t, I, k_r) qssa(R, t, I, ...
                               Ha, Hr, k_a, k_r, ...
                               alpha, gamma, delta_x);
elseif strcmp(model, "Full")
  R = [R [0 0 0 0 0 0 0 0]]; % insert starting mRNAs into molecule vector
  model = @(R, t, I, k_r) full_ode(R, t, I, ...
                                   Ha, Hr, k_a, k_r, ...
                                   k_tl, beta, P_tc, delta_x, delta_m);
else
  error("Undefined model type: %s.", model);
endif


% data scaling for easier visualization
timescale = 1e5; % 10^5 seconds
molscale = 1e-9; % nM

% run the simulation defined by the current experiment's class
if strcmp(experiment_class, "single")
  % fix period and model
  I = @(t) I(t, input_period);
  model = @(R, t) model(R, t, I, k_r);

  % run single simulatin
  step_number = ceil(simulation_total_time / simulation_step);
  Rs = euler_simulate(model, R, step_number, simulation_step);

  % get time, protein and input concentration series
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
    axis([-Inf,+Inf, 0,input_amplitude * 1.23 / molscale]);
    pbaspect([1 0.334 1]);
    legend('location', 'east');
    xlabel("Time (10^5 seconds)");
    ylabel("Concentration (nM)");
    subplot(2, 1, 1); % may or may not be needed for the model plot
  endif

  % plot model behaviour
  plot(t,R1s,'--m;R1;', t,R2s,':k;R2;', t,R3s,'-r;R3;', t,R4s,'-.g;R4;');
  axis([min(t), max(t)]);
  pbaspect([1 0.334 1]);
  legend('location', 'east');
  xlabel("Time (10^5 seconds)");
  ylabel("Concentration (nM)");

elseif strcmp(experiment_class, "period")
  % frequency response vectors
  input_periods = [];
  output_periods = [];

  T_range = period_experiment_range;
  for input_period = T_range(1) : T_range(2) : T_range(3)
    % fix this iteration's input and model
    temp_I = @(t) I(t, input_period);
    temp_model = @(R, t) model(R, t, temp_I, k_r);

    % actual simulation
    simulation_total_time = 5 * input_period;
    step_number = ceil(simulation_total_time / simulation_step);
    Rs = euler_simulate(temp_model, R, step_number, simulation_step);

    % detecting output period by via [R1] peaks
    [pks, idx] = findpeaks(Rs(:,1));
    if length(pks) < 2
      out_period = 0;
    else
      out_period = (idx(end) - idx(end-1)) * simulation_step; % last oscillation
    endif

    % store frequency response
    input_periods =  [input_periods  input_period];
    output_periods = [output_periods out_period];
  endfor

  % plot frequency response
  input_periods /= timescale;
  output_periods /= timescale;
  plot(input_periods, output_periods, '@b');
  xlabel("Input period (10^5 seconds)");
  ylabel("Output period (10^5 seconds)");

elseif strcmp(experiment_class, "oscillator")
  % frequency response vectors
  DC_range = oscillator_input_range;
  DC_range = DC_range(1) : DC_range(2) : DC_range(3);
  output_periods = [];

  for input_dc_level = DC_range
    % fix this iteration's input and model
    temp_I = @(t) input_dc_level;
    temp_model = @(R, t) model(R, t, temp_I, k_r);

    % actual simulation
    step_number = ceil(simulation_total_time / simulation_step);
    Rs = euler_simulate(temp_model, R, step_number, simulation_step);

    % detecting output period by via [R1] peaks
    [pks, idx] = findpeaks(Rs(:,1));
    if length(pks) < 2
      out_period = 0;
    else
      out_period = (idx(end) - idx(end-1)) * simulation_step; % last oscillation
    endif

    % store frequency response
    output_periods = [output_periods out_period];
  endfor

  % plot oscillatory behaviour
  DC_range /= molscale;
  output_periods /= timescale;
  plot(DC_range, output_periods, 'b');
  xlabel("Input concentration (nM)");
  ylabel("Period (10^5 seconds)");

elseif strcmp(experiment_class, "switch")
  % calculate trigger boundaries
  step_number = ceil(simulation_total_time / simulation_step);
  trigger_start = round(switch_trigger_time(1) / simulation_step);
  trigger_stop = round(switch_trigger_time(2) / simulation_step);
  prelude_steps = trigger_start;
  interlude_steps = trigger_stop - trigger_start;
  postlude_steps = step_number - (prelude_steps + interlude_steps);

  % fix input and models
  I = @(t) I(t, input_period);
  normal_model = @(R, t) model(R, t, I, k_r);
  triggered_model = @(R, t) model(R, t, I, k_r + switch_R*switch_trigger_delta);

  % 3-step simulation with different binding affinity
  prelude = euler_simulate(normal_model, R, prelude_steps, simulation_step);

  interlude = euler_simulate(triggered_model, prelude(end,:), ...
                             interlude_steps, simulation_step);

  postlude = euler_simulate(normal_model, interlude(end,:), ...
                            postlude_steps, simulation_step);

  Rs = [prelude; interlude; postlude];

  % get time, protein and input concentration series
  t = (0 : simulation_step : simulation_total_time)' / timescale;
  R1s = Rs(:,1) / molscale;
  R2s = Rs(:,2) / molscale;
  R3s = Rs(:,3) / molscale;
  R4s = Rs(:,4) / molscale;
  Is = zeros(size(t));
  for i = 1 : length(t)
    Is(i) = I(i * simulation_step) / molscale;
  endfor

  % plot model behaviour
  plot(t,R1s,'--m;R1;', t,R2s,':k;R2;', t,R3s,'-r;R3;', t,R4s,'-.g;R4;');
  xlabel("Time (10^5 seconds)");
  ylabel("Concentration (nM)");

else
  error("Undefined experiment type: %s.", experiment_class);
endif

% save plotted figure (gets put in the same folder the script runs from)
print(plot_output_filename);
