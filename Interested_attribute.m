function data_table = Interested_attribute(data_table)

% Define the column numbers to remove (e.g., columns 2 and 4)
columns_to_remove = [2,3,5,6,7,8,10,12,13,15,14,16,18,20,22,24,26,27,28,30,32,33,34,36,38,39,40,42,43,44,...
    45,46,47,48,49,50,51,...
    52,53,54,56,58,59,60,61,62,64,67,69,70,72,74,75,76,78,80,82,83,84,85,86,87,88,89,90,91,92,93,94,...
    96,97,98,99,100,101,102,103]; % Specify the column numbers to remove

% Drop the specified columns
data_table(:, columns_to_remove) = []; % Remove the specified columns