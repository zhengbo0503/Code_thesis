function validate_complex_stack()
%VALIDATE_COMPLEX_STACK - Check the complex MEX wrappers and algorithms
%
%   Usage:
%       validate_complex_stack()
%
%   Purpose:
%       VALIDATE_COMPLEX_STACK performs a short preflight check before the
%       full experiments. It calls all four complex LAPACK MEX wrappers,
%       runs one double-complex MP3JacobiSVD example and one single-complex
%       SSD example, and checks LAPACK INFO values, residuals, and
%       unitarity.
%
%   Arguments:
%       None. The four complex MEX wrappers must already be built.
%
%   Outputs:
%       None. The function raises an error if a check fails and prints a
%       completion message otherwise.
%

setup_complex_paths();
assert(ismember(version('-release'),{'2025b','2026a'}), ...
    'The complex experiments target MATLAB R2025b or R2026a (found R%s).', version('-release'));

required = {'zgesvj_mex','zgejsv_mex','cgesvj_mex','cgejsv_mex'};
for k = 1:numel(required)
    assert(exist(required{k},'file') == 3, ...
        '%s is not built. Run build_all_complex_mex first.',required{k});
end

% Double-complex wrappers and the three-precision algorithm.
A = complex_randsvd(24,12,1e6,3,'double',101);
[Uz,Sz,Vz,~,~,infoz] = zgesvj_mex(A,'G');
[Ue,Se,Ve,~,~,~,infoe] = zgejsv_mex(A);
assert(infoz == 0,'ZGESVJ returned INFO=%d.',infoz);
assert(infoe == 0,'ZGEJSV returned INFO=%d.',infoe);
check_factorisation(A,Uz,Sz,Vz,5e-10,'ZGESVJ');
check_factorisation(A,Ue,Se,Ve,5e-10,'ZGEJSV');
[U,S,V] = mposj_complex(A);
check_factorisation(A,U,S,V,5e-10,'MP3JacobiSVD');

% Single-complex wrappers and the SSD algorithm.
As = complex_randsvd(20,12,1e5,3,'single',102);
[Uc,Sc,Vc,~,~,infoc] = cgesvj_mex(As,'G');
[Uce,Sce,Vce,~,~,~,infoce] = cgejsv_mex(As);
assert(infoc == 0,'CGESVJ returned INFO=%d.',infoc);
assert(infoce == 0,'CGEJSV returned INFO=%d.',infoce);
check_factorisation(As,Uc,Sc,Vc,5e-4,'CGESVJ');
check_factorisation(As,Uce,Sce,Vce,5e-4,'CGEJSV');
[Us,Ss,Vs] = mposj_ssd_complex(As);
check_factorisation(As,Us,Ss,Vs,5e-4,'MP3JacobiSVD-SSD');

fprintf('Complex preflight validation passed under MATLAB %s.\n', ...
    version('-release'));
end

function check_factorisation(A,U,S,V,tolerance,name)
%CHECK_FACTORISATION - Check one computed SVD factorisation.
Ad = double(A);
Ud = double(U);
Sd = double(S);
Vd = double(V);
singularValues = diag(Sd);
residual = norm(Ad-Ud*Sd*Vd','fro')/norm(Ad,'fro');
unitarityU = norm(Ud'*Ud-eye(size(Ud,2)),inf);
unitarityV = norm(Vd'*Vd-eye(size(Vd,2)),inf);

assert(all(isfinite(singularValues)) && all(singularValues >= 0), ...
    '%s returned invalid singular values.',name);
assert(isfinite(residual) && residual < tolerance, ...
    '%s residual is %g.',name,residual);
assert(isfinite(unitarityU) && unitarityU < 20*tolerance, ...
    '%s U unitarity error is %g.',name,unitarityU);
assert(isfinite(unitarityV) && unitarityV < 20*tolerance, ...
    '%s V unitarity error is %g.',name,unitarityV);
end
