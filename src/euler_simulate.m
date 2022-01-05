% Simulates the model defined by function f, for m steps of dt seconds over an
% initial state X (a row vector). f(t, X) should return the gradient at time t.
%
% Returns a matrix Y, where rows correspond to protein concentrations
% at each simulation step.
function Y = euler_simulate(f, X, m, dt)

  assert(size(X, 1) == 1);
  assert(m > 0);
  assert(dt >= 1);

  n = length(X);
  Y = zeros(m, n);

  for i = 1 : m
    Y(i,:) = X;
    dX = f((i-1) * dt, X);
    X += dX * dt;
  endfor

endfunction
