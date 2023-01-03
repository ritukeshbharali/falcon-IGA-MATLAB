% BFGS solver with Wolfe line search

function globdat = bfgsSolver(props,globdat)

% Increment step and time
globdat.ts.step = globdat.ts.step + 1;

if ismember(globdat.ts.step,props.tStep.cutsteps) == 1
    globdat.ts.dt = globdat.ts.dt/props.tStep.cutsize;
end

globdat.ts.t    = globdat.ts.t + globdat.ts.dt;

fprintf('\n')
disp('========================================================')
disp(['Step ', num2str(globdat.ts.step),', Method: BFGS'])
disp(['Elements = ',num2str(globdat.mesh.numElements),...
      ', DOFs = ',num2str(globdat.ndofs)])
disp('========================================================')
disp(['Iteration    | ','Error'])

errVec = zeros(props.nlSolver.maxiter,1);

II         = speye(globdat.ndofs);

% Apply Dirichlet Boundary Condition
globdat.state = applyConstraints(globdat.state,globdat.ts.t,globdat.consTable);

% Assemble only internal force
% globdat = ModelAction(props,globdat,'AssembleIntForce');

% Set iteration count to zero and allocate vector to store error
iter   = 0;

% Enter staggered solution process
while true
    
    % Iteration counter increment, and store old internal force
        iter     = iter + 1;

        fint0    = globdat.fint;

        if iter == 1
            globdat = ModelAction(props,globdat,'AssembleDiagSystem');
        else
            %tmp       = globdat.K * s;
            %globdat.K = globdat.K + (y*y')/(y'*s) - (tmp*tmp')/(s'*tmp);
            rho   = 1.0/dot(y,s);
            tmp1  = (II-rho*s*y');
            tmp2  = rho*(s*s');
            globdat.K = tmp1*globdat.K*tmp1 + tmp2;
        end

        % fint0           = globdat.fint;

        % Current search direction
        delSol          = LinearSolver(globdat.CMat'*globdat.K*globdat.CMat, ...
                          globdat.CMat'*(-globdat.fint),props.linSolver);
        delSol          = globdat.CMat*delSol;

        % Store x0
        x0 = globdat.state;
        s0 = globdat.fint' * delSol;

        % Carry out line search
        k = 0;
        alfa = 1.0;
        alfaL = 0.1;
        alfaR = 1.0;
%         while true
%             k = k + 1;
%             globdat.state = x0 + alfa * delSol;
%             globdat = ModelAction(props,globdat,'AssembleIntForce');
%             s       = globdat.fint' * delSol;
% 
%             if abs(s) < 0.9 * abs(s0)
%                 disp(['Success! alpha = ',num2str(alfa)])
%                 break;
%             elseif k == 10
%                 disp('Failed!')
%                 alfa = 1.0;
%                 break;
%             end
% 
%             if isinf(alfaR)
%                 alfa = 2 * alfaL;
%             else
%                 alfa = 0.5 * (alfaL + alfaR);
%             end
% 
%             
%         end


        globdat.state   = x0  + alfa * delSol;

        globdat = ModelAction(props,globdat,'AssembleIntForce');
        s = delSol;
        y = globdat.fint - fint0;
                
        % Compute error and store in errCol vector
        if iter == 1
            err0 = norm(globdat.CMat'*globdat.fint);
            err  = err0;
        else
            err = norm(globdat.CMat'*globdat.fint)/max(1,err0);
        end
        errVec(iter) = err;
        
        % Display current iteration and error info
        disp([num2str(iter),'            |',num2str(err)])
        
        % Check for convergence
        if err < props.nlSolver.tol
        
            % Display message to command window            
            fprintf('\n')
            disp(['Converged in ',num2str(iter),' iteration(s).'])
            loads = abs(sum(globdat.fint(globdat.loadDofs)));
            disp(['Load = ',num2str(loads),' Displacement = ',num2str(globdat.ts.t)])
            % Store load-displacement data
            globdat.lodi(globdat.ts.step+1,:) = [loads,globdat.ts.t];
            
            %-------------------------------------------------------------%
            % Mesh refinement and solution transfer                       %
            %-------------------------------------------------------------%
            
            % Check if mesh needs to be refined
            if sum([globdat.markRef{:}]) && props.geom.AMR
                           
                % Display message to command window
                fprintf('\n')
                disp('**** MESH CHANGED, re-do step! ****')
    
                % Refine and transfer solution, build constraint matrix
                globdat = ModelAction(props,globdat,'refineNTransfer');
                                      
                % Do not accept the current load increment
                globdat.ts.step = globdat.ts.step - 1;
                globdat.ts.t    = globdat.ts.t - globdat.ts.dt;
                
                break;

            end  
            
        globdat = ModelAction(props,globdat,'commit');
        
        % Re-size error vector and store in globdat
        errVec(iter+1:props.nlSolver.maxiter) = [];
        globdat.errVec{globdat.ts.step,1} = errVec;
        
        break;
            
        %-------------------------------------------------------------%
        % Terminate NR iterations when max iter is reached or when    %
        % the error blows up                                          %
        %-------------------------------------------------------------%
        elseif iter == props.nlSolver.maxiter || err > 1000.0
            
            % Message to command window
            disp('**** WARNING: Failed to converge! ****')
            
            % Revert variables to oldStep
            globdat.state  = globdat.state0;
            
            % Reduce the stepsize for prescribed displacement method
            globdat.ts.step = globdat.ts.step - 1;
            globdat.ts.t    = globdat.ts.t - globdat.ts.t;
            globdat.ts.dt   = globdat.ts.dt/props.tStep.cutsize;
                    
            break;
            
        end

end


end % eof

