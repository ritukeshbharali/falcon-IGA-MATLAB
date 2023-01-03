function props = inputReader(casename)

% Check if the input file exists
if isfile(fullfile('./input_data/',casename,'/input_data.txt')) == 0
    error('input_data.txt file does not exist')
end

% If file exists


end