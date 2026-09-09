%TEST4 -- timing test, changing number of cols
addpath("../shared/");

close all; clear; rng(1);
fprintf("[START] svals_timing\n");

for i = 1:10
    A = gallery('randsvd', [1000,800], 1e8, 3);
    [U,S,V] = mposj(A,3);
end