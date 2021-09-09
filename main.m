clear
close
clc 

% If the data file doesn't exist then
% run the script, load the resulting data for initial values as a struct
if(isfile('data.mat'))
    inits = load('data.mat');
else
    run('givenData.m')
    inits = load('data.mat');
end
