function root = setup_complex_paths()
%SETUP_COMPLEX_PATHS - Add the complex experiment code to the MATLAB path
%
%   Usage:
%       root = setup_complex_paths()
%
%   Purpose:
%       SETUP_COMPLEX_PATHS adds the complex experiment directory, the
%       repository root containing shared routines, and the four complex
%       MEX-wrapper directories to the MATLAB search path.
%
%   Arguments:
%       None.
%
%   Outputs:
%       (1) root - Absolute path of the complex_experiments directory.
%

root = fileparts(mfilename('fullpath'));
parent = fileparts(root);

addpath(root, fullfile(parent,'shared'), ...
    fullfile(root, 'get_zgesvj'), ...
    fullfile(root, 'get_zgejsv'), ...
    fullfile(root, 'get_cgesvj'), ...
    fullfile(root, 'get_cgejsv'));
end
