%=========================================================================%
% falcon-IGA-MATLAB is an isoGeometric analysis program from the falcon
% numerical analyses suite.
% 
% Authors:     Ritukesh Bharali (Chalmers University of Technology, Sweden)
%              Somdatta Goswami (Brown University, USA)
%=========================================================================%

function globdat = main(problem)

% User prompt if no casename is provided
if nargin ~= 1
    problem = input('Enter problem to run: ','s');
end

% Clear workspace
close all;
format long;
warning('off')

% Start timer
tic

%-------------------------------------------------------------------------%
% ADD PATH TO DIRECTORIES
%-------------------------------------------------------------------------%

% Path to external libraries and source directories
addpath(genpath('./ext'));
addpath(genpath('./src/io'));
addpath(genpath('./src/constraints'));
addpath(genpath('./src/material'));
addpath(genpath('./src/solvers'));
addpath(genpath('./src/utils'));

% Remove path to femodel and input_data to avoid function clashes
rmpath(genpath('./src/fem'))
rmpath(genpath('./demo'))

% Add path to input directory
if isfolder(fullfile('./demo/',problem))
    addpath(fullfile('./demo/',problem));
else
    error(' - Problem does not exist in ./demo')
end

% Create an output directory
dir_output = sprintf('./output/%s', problem);
if ~exist(dir_output, 'dir')
    mkdir(dir_output)
end

% Remove path to output, and add only case specific directory
rmpath(genpath(dir_output))
addpath(fullfile(dir_output,problem))

%-------------------------------------------------------------------------%
% INITIALIZE
%-------------------------------------------------------------------------%

% Read input data and check its validity
props           = inputData;
props           = checkInputData(props);

% Add path to the chosen femodel
addpath(fullfile('src/fem/',props.femodel.type))

% Initialise the global database
globdat         = setupProblem(props);

%-------------------------------------------------------------------------%
% RUN
%-------------------------------------------------------------------------%

while globdat.active
           
    % Solve a step
    globdat = solveStep(props,globdat);
    
    % Check for exit
    globdat = checkExit(props,globdat);
    
    % Post-process
    globdat = postProcessStep(dir_output,props,globdat);
    
end

%-------------------------------------------------------------------------%
% SHUTDOWN
%-------------------------------------------------------------------------%

shutdownProblem(dir_output,props,globdat);
disp(' ')
disp([' - *** TIME ELAPSED = ',num2str(toc),' seconds.'])
disp(' ')
end

%=========================================================================%
% End of main function
%=========================================================================%