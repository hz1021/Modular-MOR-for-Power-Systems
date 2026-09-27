% Add all subfolders to the MATLAB path
project_folder = pwd;  % or specify your project's main folder if needed
addpath(genpath(project_folder));

% === Project Startup Script ===
% Purpose: Set up Simulink cache and codegen folders inside 'SimulinkFiles'

disp('Running startup.m: Setting Simulink cache/codegen folders...');

% Get the full path to this startup.m file
thisFilePath = mfilename('fullpath');
[projectRoot, ~, ~] = fileparts(thisFilePath);

% Define SimulinkFiles folder path
simulinkFolder = fullfile(projectRoot, 'simulink_files');

% Ensure the SimulinkFiles directory exists
if ~exist(simulinkFolder, 'dir')
    mkdir(simulinkFolder);
end

% Set simulation cache and code generation folders
Simulink.fileGenControl('set', ...
    'CacheFolder', simulinkFolder, ...
    'CodeGenFolder', simulinkFolder, ...
    'createDir', true);

disp(['Simulink cache and codegen will be stored in: ' simulinkFolder]);
