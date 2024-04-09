function samples = PoissonBetaSliceSampleBeta(logdist, xx, widths)
% SliceSampleBeta
% 
% By Jong Kyoung Kim 
% Modified from slice_sample.m by Iain Murray May 2004
% jkkim@ebi.ac.uk
% Last Update: August 31 2011

maxStepSize = 1000;

log_Px = feval(logdist, xx);

log_uprime = log(rand) + log_Px;

xprime = xx;

% Create a horizontal interval (x_l, x_r) enclosing xx
rr = rand;
x_l = xx - rr*widths;
x_r = xx + (1-rr)*widths;
J = floor(rand*maxStepSize);
K = maxStepSize-1-J;
while (feval(logdist, x_l) > log_uprime) && x_l - widths >=0 && J > 0
    x_l = x_l - widths;
    J = J - 1;
end;
while (feval(logdist, x_r) > log_uprime) && x_r + widths <=1 && K > 0
    x_r = x_r + widths;
    K = K - 1;
end;
if x_l < 0
    x_l = 0;
end;
if x_r > 1
    x_r = 1;
end;

% Propose xprimes and shrink interval until good one found
while 1
    xprime = rand*(x_r - x_l) + x_l;
    log_Px = feval(logdist, xprime);
    if log_Px > log_uprime
        break % this is the only way to leave the while loop
    else
        % Shrink in
        if xprime > xx
            x_r = xprime;
        elseif xprime < xx
            x_l = xprime;
        else
            error('BUG DETECTED: Shrunk to current position and still not acceptable.');
        end;
    end;
end;
samples = xprime;