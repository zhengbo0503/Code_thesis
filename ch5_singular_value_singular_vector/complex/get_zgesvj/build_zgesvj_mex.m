function build_zgesvj_mex()
%BUILD_ZGESVJ_MEX - Build the double-complex GESVJ MEX wrapper
%
%   Usage:
%       build_zgesvj_mex()
%
%   Purpose:
%       BUILD_ZGESVJ_MEX compiles zgesvj_mex.c with MATLAB's interleaved
%       complex API and links it against the LAPACK library supplied by
%       MATLAB R2025b. The caller's working directory is restored afterwards.
%
%   Arguments:
%       None.
%
%   Outputs:
%       None. The function writes zgesvj_mex.<mexext> to
%       complex_experiments.
%

assert(strcmp(version('-release'),'2025b'),'This build targets MATLAB R2025b.');
sourceDir = fileparts(mfilename('fullpath'));
root = fileparts(sourceDir);
startDir = pwd;
restoreDirectory = onCleanup(@() cd(startDir));
cd(root);
mex('-R2018a','-v',fullfile(sourceDir,'zgesvj_mex.c'),'-lmwlapack');

clear restoreDirectory
end
