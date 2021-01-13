% reset path
path(pathdef)

% add this utils folder
addpath(pwd)

% add SCvx files
addpath('../')

% add CVX files to path
cvx_path = 'cvx/';
addpath(cvx_path)
addpath(strcat(cvx_path,'builtins/'))
addpath(strcat(cvx_path,'commands/'))
addpath(strcat(cvx_path,'functions/'))
addpath(strcat(cvx_path,'functions/vec_'))
addpath(strcat(cvx_path,'lib/'))
addpath(strcat(cvx_path,'structures/'))
clear cvx_path

% add ECOS files to path
ecos_path = 'ecos/';
addpath(ecos_path)
addpath(strcat(ecos_path,'bin/'))
addpath(strcat(ecos_path,'src/'))
clear ecos_path