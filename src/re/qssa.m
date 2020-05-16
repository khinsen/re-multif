% Calculates the protein concentrations gradient of the network (using the QSSA
% model) for given repressor concentrations R at instant t. I(t) should evaluate
% to the input signal at time t; Ha(X) and Hr(X) must compute the Hill functions
% for activation and repression of protein X; and alpha, gamma and delta are
% reaction constants extracted from the QSSA.
function dR = qssa(R, t, I, Ha, Hr, alpha, gamma, delta)

  dR(1) = alpha*Ha(I(t))*Hr(R(2)) + gamma*Hr(R(3)) - delta*R(1);
  dR(2) = alpha*Ha(I(t))*Hr(R(4)) + gamma*Hr(R(3))*Hr(R(4)) - delta*R(2);
  dR(3) = alpha*Ha(I(t))*Hr(R(4)) + gamma*Hr(R(1)) - delta*R(3);
  dR(4) = alpha*Ha(I(t))*Hr(R(2)) + gamma*Hr(R(1))*Hr(R(2)) - delta*R(4);

endfunction
