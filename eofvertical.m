function [eof_in_vertical,pc,modevar,expvar] = eofvertical(X,delta_z)
%%% The function eofvertical computes vertical Empirical Orthogonal 
%%% Functions (EOFs) from simultaneous horizontal velocity data with 
%%% nonuniform vertical sampling. It employs the depth-averaged 
%%% velocity variance as a metric for ensuring vertical orthogonality of 
%%% the EOF modes.

%%% Input
%%% X is the timeseries of horizontal velocity at different depths where 
%%% the matrix size is [Depths × time_steps].
%%% delta_z is the layer thickness.

%%% Output
%%% eof_in_vertical returns the vertical EOFs that satisfy the orthogonality 
%%% condition based on depth-averaged velocity variance. Matrix size: [Modes × depths].
%%% pc returns the principal components. Matrix size: [Modes × time_steps]
%%% modevar returns the variance explained by each mode. Matrix size: [Modes × 1]
%%% expvar returns the percentage of total variance explained by each mode. [Modes × 1]

%%% Acknowledgement
%%% This function is inspired by the eof function in Climate Data Toolbox (CDT) developed by Chad A. Greene.
%%% The calculations were conducted with valuable guidance from my advisors, 
%%% Andrew L. Stewart and James C. McWilliams.

%%% Department of Atmospheric and Oceanic Sciences
%%% University of California, Los Angeles
%%% Cheng Yang (Sunny) Yeh
%%% July 2025

%% Check size
n = size(X);

assert(length(n)==2,'Input error: data matrix X has to be 2 dimensions');

%%% Reshape X matrix into [Depth × time_steps]
if n(2) > n(1)
    X = X;
elseif n(1) > n(2)
    X = X';
end
n = size(X);
depth_n = n(1);
time_n = n(2);

%% Calculate EOFs
W_sqrt = sqrt(1/sum(delta_z) .* diag(delta_z));
W = W_sqrt .^2;
% W_c = nan; %%% Complex eof (add this function later);

XTILDE = W_sqrt * X;
CTILDE = 1/time_n * XTILDE * XTILDE.';
[VTILDE,DTILDE] = eigs(CTILDE,depth_n);     %%% Calculate all vertial modes

%%% Recover the real-space EOFs
real_VTILDE = W^(-1/2) * VTILDE;
real_VTILDE = real_VTILDE';     %%% Matrix size: [Modes × depths]


%%% Recover the real principle component
real_PC = real_VTILDE * W * X;

%%% Output
pc = real_PC;    %%% Matrix size: [Modes × time_steps]
eof_in_vertical = real_VTILDE;  %%% Matrix size: [Modes × depths]
modevar = diag(DTILDE);     %%% Matrix size: [Modes × 1]
expvar = diag(DTILDE./sum(sum(DTILDE)));    %%% Matrix size: [Modes × 1]



end
