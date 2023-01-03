%=========================================================================%
% setupProblem (creates a globdat struct to store runtime variables)
% 
% INPUT: props - struct containing all user-defined model input parameters
%
%=========================================================================%

function globdat = setupProblem(props)

%-------------------------------------------------------------------------%
% Initialize globdat as struct
%-------------------------------------------------------------------------%

globdat = struct;


%-------------------------------------------------------------------------%
% Create mesh and compute Basis functions
%-------------------------------------------------------------------------%

if props.geom.AMR == false % (Adaptive mesh refinement switched OFF)

    load MeshInfo.mat; 
    load PHTelem.mat; 
    load controlPts.mat;

    globdat.PHTelem        = PHTelem;
    globdat.ctrlPts        = controlPts;
    globdat.mesh           = MeshInfo;

else % (Adaptive mesh refinement switched ON)

    [PHTelem,ctrlPts,mesh] = modelMesh(props.geom);

    globdat.PHTelem        = PHTelem;
    globdat.ctrlPts        = ctrlPts;
    globdat.mesh           = mesh;

end

Basis                      = cartdev(globdat.PHTelem,globdat.ctrlPts,props.geom);
globdat.Basis              = Basis;


%-------------------------------------------------------------------------%
% Initialize FE Model
%-------------------------------------------------------------------------%

globdat = ModelAction(props,globdat,'INIT');


%-------------------------------------------------------------------------%
% Initialize runtime variables (state vectors)
%-------------------------------------------------------------------------%

% Initialize state vectors (current, old step, old old step)
globdat.state        = zeros(globdat.ndofs,1);
globdat.state0       = zeros(globdat.ndofs,1);
globdat.state00      = zeros(globdat.ndofs,1);

% Initialize state increment vectors (current, old step)
globdat.Dstate       = zeros(globdat.ndofs,1);
globdat.Dstate0      = zeros(globdat.ndofs,1);

% Initialize int force vectors (current, old step)
globdat.fint         = zeros(globdat.ndofs,1);
globdat.fint0        = zeros(globdat.ndofs,1);


%-------------------------------------------------------------------------%
% Get and process model constraints
%-------------------------------------------------------------------------%

globdat = getConstraints(props,globdat);


%-------------------------------------------------------------------------%
% Initialize timestepping and solver specific runtime variables
%-------------------------------------------------------------------------%

% Parameters for timestepping
globdat.ts.step         = 0;
globdat.ts.t            = 0;
globdat.ts.dt           = props.tStep.initstepsize;

% Previous stepsizes (old, oldold) for quasi-nonlinear solver
if strcmp(props.nlSolver.type,'qNonlin')
    globdat.ts.dt0      = props.tStep.initstepsize;
    globdat.ts.dt00     = props.tStep.initstepsize;
    globdat.Exstate     = zeros(globdat.ndofs,1);
end

% Some variables required for arc-length solver
if strcmp(props.nlSolver.type,'nonlinArcLength')
    globdat.method      = 'dispCtrl';
    globdat.ts.t0       = 0.0;
    globdat.ts.fail     = 0;
end

% Line-search parameter
globdat.beta            = 1.0;


%-------------------------------------------------------------------------%
% Initialize additional runtime variables for post-processing
%-------------------------------------------------------------------------%

% Array for load-displacement data
globdat.lodi            = zeros(props.tStep.nsteps+1,2);

% Create a cell to store error vectors for each step
globdat.errVec = cell(props.tStep.nsteps,1);


%-------------------------------------------------------------------------%
% Activate globdat to start timestepping
%-------------------------------------------------------------------------%

globdat.active          = true;

end
%=========================================================================%
% End of function
%=========================================================================%
