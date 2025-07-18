function recon_X = reofvertical(eof_in_vertical,pc,modes)
%%% This function reconstructs the flow using principal components,
%%% vertical EOF structures, and the specified number of modes.

%%% Input
%%% eof_in_vertical — vertical EOFs structure; matrix size: [Modes × depths]
%%% pc              — principal components; matrix size: [Modes × time_steps]
%%% modes           — mode index or range (e.g., 1 or 1:2 or 1:length(eof_in_vertical))

%%% Output
%%% recon_X, Reconstructed flow X; matrix size: [Depths × time_steps]

%%% Acknowledgement
%%% This function is inspired by the eof function in Climate Data Toolbox (CDT) developed by Chad A. Greene.
%%% The calculations were conducted with valuable guidance from my advisors, 
%%% Andrew L. Stewart and James C. McWilliams.

%%% Department of Atmospheric and Oceanic Sciences
%%% University of California, Los Angeles
%%% Cheng Yang (Sunny) Yeh
%%% July 2025

%% Check size
n = size(eof_in_vertical);
assert(length(n)==2,'Input error: eof_in_vertical has to be 2 dimensions');
assert(length(modes)<=length(eof_in_vertical) && max(modes)<=length(eof_in_vertical),'Input error: mode(s) has(have) to be no larger than the number of available modes.')

%%% Reshape pc into [Modes * time_steps]
if size(pc,2) > size(pc,1)
    pc = pc;
elseif size(pc,1) > size(pc,2)
    pc = pc';
end
n_pc = size(pc);
depth_n_pc = n_pc(1);
time_n_pc = n_pc(2);

%% Reconstruct the flow
recon_X = nan(n_pc);
for loop_X_depths = 1 : depth_n_pc
    temp = zeros(1,length(recon_X));
    for loop_modes = modes
        temp = temp + eof_in_vertical(loop_modes,loop_X_depths).*pc(loop_modes,:);  %%% eof_in_vertical size has to be [Modes * depths]
    end
    recon_X(loop_X_depths,:) = temp;
end
end
