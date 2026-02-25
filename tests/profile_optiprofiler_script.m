clear all
options.dim = 'small';
options.plibs = 's2mpj';
options.feature_name = 'plain';
% options.n_jobs = 1;

% options.solver_names = {'cbds', 'ds'};

% profile_optiprofiler(options);

% options.solver_names = {'cbds', 'newuoa'};

% profile_optiprofiler(options);

% options.solver_names = {'cbds', 'fd-bfgs'};

% profile_optiprofiler(options);

options.solver_names = {'cbds', 'cbds-orig-direction-set-from-x0'};
profile_optiprofiler(options);