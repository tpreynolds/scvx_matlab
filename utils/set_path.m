% reset path
path(pathdef)

% add this utils folder
addpath(pwd)

% add SCvx files
addpath('../')

% add CVX files to path
cvx_path = '../../../matlab/tools/CVX/';
addpath(cvx_path)
addpath(strcat(cvx_path,'builtins/'))
addpath(strcat(cvx_path,'commands/'))
addpath(strcat(cvx_path,'functions/'))
addpath(strcat(cvx_path,'functions/vec_'))
addpath(strcat(cvx_path,'lib/'))
addpath(strcat(cvx_path,'structures/'))
clear cvx_path