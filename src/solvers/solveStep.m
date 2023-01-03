%=========================================================================%
% solveStep (solves the linear/nonlinear problem with user-specified solver
%            type)
% 
% Available solvers: 
%                  * staggered     (staggered solver)
%                  * acclStaggered (Anderson accelerated staggered solver)
%                  * nonlin        (Conventional Newton-Raphson solver)
%                  * qNonlin       (quasi-NR solver)
%                  * nonlinArcLength (NR solver with arc-length control)
%=========================================================================%

function globdat = solveStep(props,globdat)

switch props.nlSolver.type
    
    case 'staggered'
        
        globdat = staggeredSolver(props,globdat);

    case 'acclStaggered'
        
        globdat = acclStaggeredSolver(props,globdat);
        
    case 'nonlin'
        
        globdat = nonlinSolver(props,globdat); 
        
    case 'qNonlin'
        
        globdat = qNonlinSolver(props,globdat); 

    case 'inflatedNonlin'

        globdat = inflatedNonlinSolver(props,globdat); 

    case 'bfgs'
        
        globdat = bfgsSolver(props,globdat); 

    case 'nonlinArcLength'

        globdat = nonlinArcLength(props,globdat);
        
    otherwise
        
        error('Wrong option!')
        
end

end

