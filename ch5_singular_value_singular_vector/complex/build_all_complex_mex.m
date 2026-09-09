function build_all_complex_mex()
%BUILD_ALL_COMPLEX_MEX - Build all four complex LAPACK MEX wrappers
%
%   Usage:
%       build_all_complex_mex()
%
%   Purpose:
%       BUILD_ALL_COMPLEX_MEX builds the ZGESVJ, ZGEJSV, CGESVJ, and CGEJSV
%       wrappers against the LAPACK library supplied by MATLAB R2025b. The
%       caller's working directory and MATLAB path are restored afterwards.
%
%   Arguments:
%       None.
%
%   Outputs:
%       None. The four MEX binaries are written to complex_experiments.
%

originalDir = pwd;
restoreOriginalDir = onCleanup(@() cd(originalDir));
root = fileparts(mfilename('fullpath'));
builderDirectories = {
    fullfile(root,'get_zgesvj')
    fullfile(root,'get_zgejsv')
    fullfile(root,'get_cgesvj')
    fullfile(root,'get_cgejsv')
    };
originalPath = path;
restorePath = onCleanup(@() path(originalPath));
addpath(builderDirectories{:});
builders = {@build_zgesvj_mex,@build_zgejsv_mex, ...
    @build_cgesvj_mex,@build_cgejsv_mex};

for k = 1:numel(builders)
    builders{k}();
    cd(originalDir);
end

clear restorePath restoreOriginalDir
end
