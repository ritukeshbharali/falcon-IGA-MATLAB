%=========================================================================%
% ModelCheckInput (performs a check on the user-supplied input parameters)
% 
% INPUT: props - struct containing all user-defined model input parameters
%
%=========================================================================%

function props = checkInputData(props)

disp(' - Model Summary')

%-------------------------------------------------------------------------%
% Check if FEModel type exists in the software library
%-------------------------------------------------------------------------%

if not(isfolder(fullfile('./src/fem/',props.femodel.type)))
    fullfile('./src/fem/',props.femodel.type)
    error('FE Model does not exist! Available: Elasticity, PhaseFieldFracture')
else
    disp(['   + FE Model: ',props.femodel.type])
end

%-------------------------------------------------------------------------%
% Check if the minimal material properties are present to run the model
%-------------------------------------------------------------------------%

% Check if viscous regularization is present in Phase-field models
if strcmp(props.femodel.type,'PhaseField') == 1 ...
          && isfield(props.mat,'visc') == 0
    props.mat.visc = 0;
    disp('   + Viscous regularization is absent!')
end

if strcmp(props.femodel.type,'AmorPhaseField') == 1 ...
          && isfield(props.mat,'visc') == 0
    props.mat.visc = 0;
    disp('   + Viscous regularization is absent!')
end



end