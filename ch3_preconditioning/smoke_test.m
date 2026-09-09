function nfail = smoke_test()
%SMOKE_TEST - Chapter 3 smoke test
%
%   Usage:
%       nfail = smoke_test()
%
%   Purpose:
%       Exercises every function in ch3_preconditioning/shared on small
%       inputs and checks the defining property of each. Does not run the
%       full experiments. Returns the number of failed checks.
%
%       Advanpix (mp) is optional here: checks that need it are skipped
%       with a notice if it is not on the path.
%
%   Author:
%       Smoke test for the thesis code reorganisation.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'shared'));
R = {};

fprintf('=== Chapter 3: preconditioning ===\n');
rng(0);
n = 40;
A = gallery('randsvd', n, -1e6, 3);   % SPD
G = randn(60, 25);                    % general, tall

% --- float_params -------------------------------------------------------
try
    [ud,~,~,~,pd] = float_params('d');
    [us,~,~,~,ps] = float_params('s');
    uq            = float_params('q');
    R = chk(R,'float_params d/s/q ordering', uq < ud && ud < us, ...
        sprintf('uq=%.3e ud=%.3e us=%.3e', uq, ud, us));
    R = chk(R,'float_params double u = 2^-53', abs(ud - 2^-53) < eps, sprintf('%.17g', ud));
    R = chk(R,'float_params precisions', pd == 53 && ps == 24, sprintf('pd=%d ps=%d', pd, ps));
catch e
    R = bad(R,'float_params', e);
end

% --- off ----------------------------------------------------------------
try
    o1 = off(diag(rand(n,1)));
    o2 = off(A);
    R = chk(R,'off(diagonal) == 0', o1 < 1e-14, sprintf('%.3e', o1));
    R = chk(R,'off(full) > 0', o2 > 0, sprintf('%.3e', o2));
catch e
    R = bad(R,'off', e);
end

% --- orthogonalisers: mgs, mgs2, mgs3, ns ----------------------------
try
    [Q,Rr] = mgs(G);
    R = chk(R,'mgs orthogonality', norm(Q'*Q - eye(size(Q,2)),'fro') < 1e-10, ...
        sprintf('%.3e', norm(Q'*Q - eye(size(Q,2)),'fro')));
    R = chk(R,'mgs reconstruction', relerr(G, Q*Rr) < 1e-12, sprintf('%.3e', relerr(G,Q*Rr)));
catch e
    R = bad(R,'mgs', e);
end
try
    [Qt,Rr2] = mgs3(G, 5);
    R = chk(R,'mgs3 orthogonality', norm(Qt'*Qt - eye(size(Qt,2)),'fro') < 1e-10, ...
        sprintf('%.3e', norm(Qt'*Qt - eye(size(Qt,2)),'fro')));
    R = chk(R,'mgs3 reconstruction', relerr(G, Qt*Rr2) < 1e-12, sprintf('%.3e', relerr(G,Qt*Rr2)));
catch e
    R = bad(R,'mgs3', e);
end
try
    [U,~,k] = ns(orthish(n));
    R = chk(R,'ns orthogonality', norm(U'*U - eye(n),'fro') < 1e-8, ...
        sprintf('%.3e after %d its', norm(U'*U-eye(n),'fro'), k));
catch e
    R = bad(R,'ns', e);
end
% --- scond: the merged function, both scalings --------------------------
try
    sD = scond(A,'D');
    sC = scond(G,'C');
    sAuto = scond(A);          % symmetric -> should pick 'D'
    sGauto = scond(G);         % nonsymmetric -> should pick 'C'
    R = chk(R,'scond D finite/positive', isfinite(sD) && sD > 0, sprintf('%.3e', sD));
    R = chk(R,'scond C finite/positive', isfinite(sC) && sC > 0, sprintf('%.3e', sC));
    R = chk(R,'scond auto == D for symmetric',  abs(sAuto - sD)/sD < 1e-12, sprintf('%.3e vs %.3e', sAuto, sD));
    R = chk(R,'scond auto == C for general',    abs(sGauto - sC)/sC < 1e-12, sprintf('%.3e vs %.3e', sGauto, sC));
    R = chk(R,'scond(A,D) <= cond(A) for SPD',  sD <= cond(A)*(1+1e-8), sprintf('scond=%.3e cond=%.3e', sD, cond(A)));
    ok = false;
    try, scond(A,'X'); catch, ok = true; end
    R = chk(R,'scond rejects bad scaling type', ok, '');
catch e
    R = bad(R,'scond', e);
end

% --- scaled_cond --------------------------------------------------------
try
    r = scaled_cond(A);
    R = chk(R,'scaled_cond finite/positive', isfinite(r) && r > 0, sprintf('%.3e', r));
catch e
    R = bad(R,'scaled_cond', e);
end

% --- mpdiag (needs Advanpix) -------------------------------------------
if have('mp')
    try
        [Qt2, Dt] = mpdiag(A);
        R = chk(R,'mpdiag orthogonality', norm(double(Qt2)'*double(Qt2) - eye(n),'fro') < 1e-8, ...
            sprintf('%.3e', norm(double(Qt2)'*double(Qt2)-eye(n),'fro')));
        R = chk(R,'mpdiag returns diagonal', isvector(Dt) || isdiag(double(Dt)), '');
    catch e
        R = bad(R,'mpdiag', e);
    end
else
    fprintf('  SKIP mpdiag (Advanpix mp not on path)\n');
end

% --- mywritetable -------------------------------------------------------
try
    tmp = [tempname '.csv'];
    mywritetable(rand(4,3), {'a','b','c'}, tmp);
    T = readtable(tmp);
    R = chk(R,'mywritetable roundtrip', height(T)==4 && width(T)==3 && ...
        isequal(T.Properties.VariableNames, {'a','b','c'}), sprintf('%dx%d', height(T), width(T)));
    delete(tmp);
catch e
    R = bad(R,'mywritetable', e);
end

nfail = summary(R, 'Chapter 3');
end

% ---------------------------------------------------------------- helpers
function R = chk(R, name, cond, detail)
R{end+1} = struct('name',name,'pass',logical(cond),'detail',detail);
fprintf('  [%s] %-42s %s\n', tf(cond), name, detail);
end
function R = bad(R, name, e)
R{end+1} = struct('name',name,'pass',false,'detail',e.message);
fprintf('  [FAIL] %-42s ERROR: %s\n', name, e.message);
end
function s = tf(c), if c, s='PASS'; else, s='FAIL'; end, end
function r = relerr(A,B), r = norm(A-B,'fro')/max(norm(A,'fro'),realmin); end
function t = have(n), t = ~isempty(which(n)); end
function Q = orthish(n)
[Q,~] = qr(randn(n)); Q = Q + 1e-3*randn(n)/n;   % near-orthogonal
end
function nfail = summary(R, label)
p = sum(cellfun(@(x)x.pass, R)); t = numel(R); nfail = t - p;
fprintf('  ---- %s: %d/%d passed, %d failed ----\n', label, p, t, nfail);
end
