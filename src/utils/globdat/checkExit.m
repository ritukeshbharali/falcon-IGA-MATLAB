%=========================================================================%
% checkExit (sets the globdat.active to false if the time-stepping needs to
%            be terminated)
% 
% Available criterion for termination: 
%                  * total number of steps reached
%                  * current load / peak load < user-defined number
%=========================================================================%

function globdat = checkExit(props,globdat)

if globdat.ts.step == props.tStep.nsteps || ...
   globdat.lodi(globdat.ts.step+1,1)/max(globdat.lodi(:,1)) < ...
   props.tStep.finalLoadRatio

% Set globdat.active to false to exit time-stepping
globdat.active = false;

end

end

