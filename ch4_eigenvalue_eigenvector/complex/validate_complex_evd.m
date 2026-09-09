function report = validate_complex_evd(outputRoot)
%VALIDATE_COMPLEX_EVD - Small numerical checks without running the full grids
%   Usage: report = validate_complex_evd(outputRoot)
%   Purpose: Check rotation phases, real reduction, residuals, unitarity,
%       scaled inputs and both reduced experiment/export paths.
%   Arguments: outputRoot defaults to a new temporary directory.
%   Outputs: report contains check counts and largest residual/unitarity.
root = setup_complex_evd_paths();
if nargin < 1, outputRoot = tempname; end
if ~isfolder(outputRoot), mkdir(outputRoot); end
assert(exist('mp','class') == 8,'Add the Advanpix toolbox to the path.');
checks = 0; maxResidual = 0; maxUnitarity = 0;
endpoint = repmat(struct('mode',NaN,'forward',NaN,'residual',NaN, ...
    'unitarity',NaN,'scaledCondition',NaN,'realCaseReference',NaN),1,5);
digitsBefore = mp.Digits();

% A phase error in J can pass real-input checks: exercise arbitrary phases.
for diagonal = [1 2 3]
    for phase = [1 -1 1i -1i exp(0.37i) exp(-2.1i)]
        A = [diagonal,0.4*phase;0.4*conj(phase),2];
        [c,s] = jacobi_pair_complex(A,1,2);
        J = [c,s;-conj(s),c]; B = J'*A*J;
        assert(norm(J'*J-eye(2),'fro') < 20*eps);
        assert(abs(B(1,2)) < 20*eps*norm(A,'fro'));
        [V,D,~,~,info] = cjacobi_complex(A);
        [f,r,o] = compute_error_evals_complex(A,V,D,reference_eigenvalues_complex(A));
        assert(info == 0 && f < 1e-14 && r < 1e-14 && o < 1e-14);
        maxResidual = max(maxResidual,r); maxUnitarity = max(maxUnitarity,o);
        checks = checks+1;
    end
end

% Zero off-diagonals, n=1, and eigenvalues small enough to expose rounding.
for diagonal = {2,[1;2^-40;2^-100]}
    A = diag(diagonal{1});
    [~,D,nrot,~,info] = cjacobi_complex(A);
    assert(info == 0 && nrot == 0 && isequal(D,sort(diagonal{1},'descend')));
    checks = checks+1;
end

% All five spectra; reproducibility uses a private stream, preserving rng.
state = rng;
for mode = 1:5
    for kappa = [1e2 1e6 1e10 1e15]
        [A,~,metadata] = complex_randhpd(12,kappa,mode,'double',100+mode);
        assert(ishermitian(A) && ~isreal(A));
        assert(isequal(A,complex_randhpd(12,kappa,mode,'double',100+mode)));
        assert(isequal(state,rng) && metadata.prescribed_kappa <= kappa*(1+eps));
        dref = reference_eigenvalues_complex(A);
        for method = {'jacobi','mp2','mp3'}
            if strcmp(method{1},'jacobi')
                [V,D,~,~,info] = cjacobi_complex(A);
            else
                [V,D,~,~,reference,scaledCondition,info] = mp_pjacobi_complex(A,string(method{1}));
            end
            [f,r,o] = compute_error_evals_complex(A,V,D,dref);
            assert(info == 0 && all(isfinite([f r o])));
            assert(r < 2e-13 && o < 2e-13,'Residual or unitarity check failed.');
            if strcmp(method{1},'mp3')
                if kappa <= 1e10
                    assert(f < 1e-12,'MP3 moderate-condition accuracy check failed.');
                else
                    % At the endpoint, the preconditioned scaled condition can
                    % be large. Record accuracy instead of assuming a fixed floor.
                    endpoint(mode) = struct('mode',mode,'forward',f, ...
                        'residual',r,'unitarity',o, ...
                        'scaledCondition',double(scaledCondition), ...
                        'realCaseReference',double(reference));
                    fprintf('ENDPOINT mode=%d forward=%.3e scaledCondition=%.3e\n', ...
                        mode,f,double(scaledCondition));
                end
            end
            maxResidual = max(maxResidual,r); maxUnitarity = max(maxUnitarity,o);
            checks = checks+1;
        end
    end
end

% Agreement with the real source functions on well-conditioned real inputs.
stream = RandStream('mt19937ar','Seed',42); X = randn(stream,10);
A = X'*X+eye(10);
for method = {'jacobi','mp2','mp3'}
    if strcmp(method{1},'jacobi')
        [Vr,Dr] = cjacobi(A); [Vc,Dc] = cjacobi_complex(A);
    else
        [Vr,Dr] = mp_pjacobi(A,string(method{1}));
        [Vc,Dc] = mp_pjacobi_complex(A,string(method{1}));
    end
    assert(norm(Dr-Dc,inf)/norm(Dr,inf) < 2e-14);
    assert(norm(A-Vr*diag(Dr)*Vr','fro')/norm(A,'fro') < 2e-13);
    assert(norm(A-Vc*diag(Dc)*Vc','fro')/norm(A,'fro') < 2e-13);
    checks = checks+1;
end

% Exercise the original MP scaling paths, with no change to their formulas.
A = complex_randhpd(6,1e3,3,'double',4);
for factor = [1e-200 1e200]
    B = factor*A;
    [V,D,~,~,~,~,info] = mp_pjacobi_complex(B,"mp3");
    dref = reference_eigenvalues_complex(B);
    f = double(max(abs(dref-D)./abs(dref)));
    r = norm(A-V*diag(double(D)/factor)*V','fro')/norm(A,'fro');
    assert(info == 0 && f < 1e-12 && r < 2e-13);
    checks = checks+1;
end

% Both reduced sweeps execute the same method/export code as the full grid.
resultsK = test_complex_evals_varying_kappa('smoke',outputRoot);
resultsN = test_complex_evals_varying_n('smoke',outputRoot);
for pair = { {resultsK,'varying_kappa'}, {resultsN,'varying_n'} }
    results = pair{1}{1}; kind = pair{1}{2};
    assert(results.passed && all(results.info(:) == 0));
    csv = fullfile(outputRoot,'data',['complex_evals_' kind '_smoke.csv']);
    values = readmatrix(csv);
    assert(isequal(size(values),[10 29]) && all(isfinite(values(:))));
    assert(isequal(values(:,6),results.actualKappas(:)));
    assert(isequal(values(:,7:10),reshape(results.forward,[],4)));
    assert(isequal(values(:,11:14),reshape(results.residual,[],4)));
    assert(isequal(values(:,15:18),reshape(results.unitarityV,[],4)));
    assert(isequal(values(:,19:21),reshape(results.info,[],3)));
    assert(isequal(values(:,22:24),reshape(results.sweeps,[],3)));
    assert(isequal(values(:,25:27),reshape(results.rotations,[],3)));
    assert(isequal(values(:,28:29),[results.scalecondAt(:),results.boundAt(:)]));
    saved = load(fullfile(outputRoot,'data',['complex_evals_' kind '_smoke.mat']));
    assert(isequaln(saved.results,results));
    for ext = {'.png','.pdf'}
        assert(isfile(fullfile(outputRoot,'plots',['complex_evals_' kind '_smoke' ext{1}])));
    end
    checks = checks+1;
end
assert(mp.Digits() == digitsBefore,'Global multiprecision setting changed.');
report = struct('passed',true,'checks',checks,'maxResidual',maxResidual, ...
    'maxUnitarity',maxUnitarity,'endpoint',endpoint,'outputRoot',outputRoot,'codeRoot',root);
save(fullfile(outputRoot,'validation.mat'),'report');
fprintf('COMPLEX_EVD_VALIDATION_PASSED checks=%d maxResidual=%.3e maxUnitarity=%.3e\n', ...
    checks,maxResidual,maxUnitarity);
end
