%=========================================================================%
% Input file for Eccentric Hole specimen (EH)
%=========================================================================%

function props = inputData()

% Initialize props as a struct
props = struct;

%-------------------------------------------------------------------------%
% FE Model
%-------------------------------------------------------------------------%

props.femodel.type = 'PhaseFieldFracture';
props.femodel.dofs = {'dx','dy','phi'};

%-------------------------------------------------------------------------%
% Geometric parameters
%-------------------------------------------------------------------------%

% Use C1 cubics with adaptivity
geom = struct;
geom.dim             = 2;
geom.nstress         = 3;
geom.patchBoundaries = {1,4,2,2; 4,3,4,4;[3,1],2,[2,4],[4,2]};

% Consider cubic B-Splines
geom.p           = 3; 
geom.q           = 3; 
geom.L           = 1; 
geom.W           = 1; 
geom.numPatches  = 4; 
geom.numElemU    = 3; 
geom.numElemV    = 3;
geom.toler       = 1e-6;
geom.ngaussX     = geom.p+1;
geom.ngaussY     = geom.q+1;
geom.maxRefLevel = 4; 
geom.threshPhi   = 0.4; 
geom.AMR         = true;

props.geom       = geom;


%-------------------------------------------------------------------------%
% Material parameters
%-------------------------------------------------------------------------%

% Elasticity parameters

plane_stress = false;
E            = 210e3;
nu           = 0.3;
lambda       = E*nu/((1+nu)*(1-2*nu));
mu           = E/(2*(1+nu));

if plane_stress == true
    lambda = (2*lambda*mu)/(lambda+2*mu);
end

mat        = struct;
mat.E      = E;
mat.nu     = nu;
mat.lambda = lambda;
mat.mu     = mu;
mat.K      = lambda + 2/3*mu;

mat.C      = zeros(3,3);
mat.C(1,1) = lambda + 2*mu;
mat.C(2,2) = lambda + 2*mu;
mat.C(3,3) = mu;
mat.C(1,2) = lambda;
mat.C(2,1) = lambda;

% Fracture parameters

mat.gc     = 2.7; %2.7; 
mat.l0     = 0.02; 
mat.eSplit = 'NoSplit';
mat.hybrid = 'false';
mat.pen    = 0.0;

% AT2 model
mat.p      = 2.0;
mat.a1     = 2.0;
mat.a2     = -0.5;
mat.a3     = 0.0;
mat.calpha = 2;
mat.eta    = 0;
mat.Psi0   = 0;

mat.visc = 0.0;

props.mat   = mat;


%-------------------------------------------------------------------------%
% Constraints
%-------------------------------------------------------------------------%

% Should be struct format

cons(1).label = {'Bottom'};
cons(1).dofs  = {'dx','dy'};
cons(1).vals  = {0,0};

cons(2).label = {'Top'};
cons(2).dofs  = {'dx','dy'};
cons(2).vals  = {0,1};

props.cons    = cons;


%-------------------------------------------------------------------------%
% Time-stepping parameters
%-------------------------------------------------------------------------%

tStep.nsteps         = 900; 
tStep.initstepsize   = 1e-4;
tStep.finalstepsize  = 1e-8;
tStep.cutsteps       = [55]; %[56];
tStep.cutsize        = 50;
tStep.finalLoadRatio = 0.02;

props.tStep          = tStep;


%-------------------------------------------------------------------------%
% Solver parameters
%-------------------------------------------------------------------------%

nlSolver.type          = 'qNonlin';%'nonlinArcLength';
nlSolver.tol           = 1e-3;
nlSolver.maxiter       = 1200;

nlSolver.arcLengthBeta       = 1.0;
nlSolver.arcLengthSwitchE    = 1e-4;
nlSolver.arcLengthMaxdTau    = 0.0125;
nlSolver.arcLengthOptIter    = 25;

linSolver.type         = 'direct';
linSolver.tol          = []; % Tolerance for iterative solvers
linSolver.prec         = []; % Preconditioner type for iterative solvers
linSolver.restart      = []; % Restart for iterative solvers
linSolver.maxiter      = []; % Maximum iterations for iterative solvers

props.nlSolver         = nlSolver;
props.linSolver        = linSolver;

%-------------------------------------------------------------------------%
% Post-processing parameters
%-------------------------------------------------------------------------%

postProc.lodi.xmax       = 0.008;
postProc.lodi.ymax       = 1000;
postProc.movie.framerate = 5; 
postProc.onScreen        = inf;
postProc.printVTK        = 5;

props.postProc           = postProc;

%=========================================================================%
% End of input file
%=========================================================================%

end