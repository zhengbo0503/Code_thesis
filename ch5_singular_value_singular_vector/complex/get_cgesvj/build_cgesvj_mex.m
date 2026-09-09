function build_cgesvj_mex()
%BUILD_CGESVJ_MEX - Build the single-complex GESVJ MEX wrapper
%
%   Usage:
%       build_cgesvj_mex()
%
%   Purpose:
%       BUILD_CGESVJ_MEX compiles cgesvj_mex.c with MATLAB's interleaved
%       complex API and links it against the LAPACK library supplied by
%       MATLAB R2025b. The caller's working directory is restored afterwards.
%
%   Arguments:
%       None.
%
%   Outputs:
%       None. The function writes cgesvj_mex.<mexext> to
%       complex_experiments.
%

assert(strcmp(version('-release'),'2025b'),'This build targets MATLAB R2025b.');
sourceDir = fileparts(mfilename('fullpath'));
root = fileparts(sourceDir);
startDir = pwd;
restoreDirectory = onCleanup(@() cd(startDir));
cd(root);
mex('-R2018a','-v',fullfile(sourceDir,'cgesvj_mex.c'),'-lmwlapack');

clear restoreDirectory
end
