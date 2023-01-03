%=========================================================================%
% shutdownProblem (carries out some post-processing steps after
%                  globdat.active is set to false, i.e., time-stepping
%                  process is terminated)
% 
% Available options: 
%                  * re-size pre-allocated load-displacement vector and
%                    print to txt file
%                  * print stepwise convergence data to files
%                  * zip the entire output folder
%=========================================================================%

function shutdownProblem(dir_output,props,globdat)

% Re-size lodi array and export to txt file
globdat.lodi(globdat.ts.step+2:end,:) = [];
writematrix(globdat.lodi,fullfile(dir_output,'lodi.txt'))  

% Export error vectors to txt file
for i = 1:length(globdat.errVec)
    writematrix(globdat.errVec{i},fullfile(dir_output,sprintf('errVec%g', i)))  
end

% Save input parameters to mat file
save(fullfile(dir_output,'inputData.mat'), 'props');

% Zip up the output folder
% zip([dir_output,'.zip'],dir_output)

end