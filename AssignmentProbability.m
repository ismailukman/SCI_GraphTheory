% % Define the size of the array
% nodeCount = 200;
% communityCount = 76;
% totalIterations = 100;
% 
% % Assuming you have the communityAssignments array already available
% 
% % Initialize variables S_g1_temp(:,:,1)
% sameCommunityCount = 0;
% 
% for k = 1:totalIterations
% % Iterate over each node pair (i, j)
%     for i = 1:nodeCount
%         for j = i+1:communityCount % Avoid duplicate pairs
%         % Check if nodes i and j are assigned to the same community in each iteration
%             if S_g1_temp(i, :, k) == S_g1_temp(j, :, k) 
%             % if S_g1_temp(5, 5, k) == S_g1_temp(5, 5, k+1) 
%                 fprintf('%d == %d ?\n', S_g1_temp(i, :, k), S_g1_temp(j, :, k) );
%                 sameCommunityCount = sameCommunityCount + 1;
%             end
%         end
%     end
% end
% 
% % Calculate probability
% totalPairs = 7 * (7 - 1) / 2; % Total number of unique node pairs communitycount*(cc-1)/2
% probability = sameCommunityCount / (totalPairs * totalIterations);
% 
% % Display result
% fprintf('Probability of node i and j being in the same community: %.4f\n', probability);

% *************************************************


% Set the number of Monte Carlo iterations
numIterations = 50;

% Initialize p-value difference matrix
p_value_diff = zeros(200, 200); 

% Generate some sample data for demonstration
% labels_s1 = randn(200, 76, 50); % Sample data for labels_s1
% labels_s2 = randn(200, 61, 50); % Sample data for labels_s2

labels_s1 = S_g1_temp;
labels_s2 = S_g2_temp; 

% Generate some sample data for demonstration
% labels_s1 = randi([0, 1], 200, 76, 50);  % Sample data for labels_s1 (binary labels)
% labels_s2 = randi([0, 1], 200, 61, 50);  % Sample data for labels_s2 (binary labels)

% Loop through each pair of variables (i, j)
for i = 1:200
    for j = i+1:200
        % Perform statistical test (e.g., t-test) between labels for variable i and j
        p_values = zeros(numIterations, 1);
        for iter = 1:numIterations
            % Randomly select indices for each label set
            idx1 = randperm(76);
            idx2 = randperm(61);
            
            % Extract labels for variable i and j from each label set
            labels_i_s1 = labels_s1(i, idx1(1:50), iter);
            labels_j_s1 = labels_s1(j, idx1(1:50), iter);
            labels_i_s2 = labels_s2(i, idx2(1:50), iter);
            labels_j_s2 = labels_s2(j, idx2(1:50), iter);
            
            % Perform the statistical test (e.g., t-test)
            [~, p_values(iter)] = ttest2(labels_i_s1, labels_j_s1); % Compare labels_s1
            [~, p_values(iter)] = ttest2(labels_i_s2, labels_j_s2); % Compare labels_s2
        end
        
        % Calculate p-value difference
        p_value_diff(i, j) = mean(abs(diff(p_values)));
        p_value_diff(j, i) = p_value_diff(i, j); % Symmetric matrix
    end
end

% Display the results
disp('P-value difference matrix:');
disp(p_value_diff);


cd '/Users/ismaila/Documents/C-Codes/SCI_GraphTheory/sci_data/SCI/modularity_var/';

filename = sprintf('p_value_diff.mat'); save(filename,'p_value_diff');
