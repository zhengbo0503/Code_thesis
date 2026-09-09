function build_cgejsv_mex()
%BUILD_CGEJSV_MEX - Build the single-complex GEJSV MEX wrapper
%
%   Usage:
%       build_cgejsv_mex()
%
%   Purpose:
%       BUILD_CGEJSV_MEX compiles cgejsv_mex.c with MATLAB's interleaved
%       complex API and links it against the LAPACK library supplied by
%       MATLAB R2025b. The caller's working directory is restored afterwards.
%
%   Arguments:
%       None.
%
%   Outputs:
%       None. The function writes cgejsv_mex.<mexext> to
%       complex_experiments.
%

assert(strcmp(version('-release'),'2025b'),'This build targets MATLAB R2025b.');
sourceDir = fileparts(mfilename('fullpath'));
root = fileparts(sourceDir);
startDir = pwd;
restoreDirectory = onCleanup(@() cd(startDir));
cd(root);
mex('-R2018a','-v',fullfile(sourceDir,'cgejsv_mex.c'),'-lmwlapack');

clear restoreDirectory
end
