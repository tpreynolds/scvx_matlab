% start fresh
clear; close all

% set path
run('../utils/set_path.m')

% startup CVX
cvx_startup

% build ecos if needed
% if ecos.o doesn't exist, run makemex there to build it