function [sys_data, branch_data] = loaddata(filename)
warning off
% Read the data from the specified sheet into a table
sys_data = table2array(readtable(filename, 'Sheet', 'system_data'));
branch_data = table2array(readtable(filename, 'Sheet', 'branch_data'));

end