function data=filter_spikes(data0)

% Define the OD threshold for change over 5 min period
threshold = 0.02;

% Length of data vector
n = length(data0);

% Iterate over the data vector to find and adjust for large changes
for i = 1:n-1
    if abs(data0(i+1) - data0(i)) > threshold
        % Calculate the shift
        shift = data0(i+1) - data0(i);
        
        % Adjust the rest of the vector to disconsider the shift
        data0(i+1:end) = data0(i+1:end) - shift;
    end
end
data=data0;
end