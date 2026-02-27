clc; clear; close all;

%% TEST CASE: 3D "Deep Trap" Basin
% -------------------------------------------------------------------------
% Objective:
%   To demonstrate failure of Default Strategy in a landscape where local 
%   gradients OVERWHELM global gradients at the micro-scale.
%
% Configuration:
%   - Global Gradient: ~ 20,000 (Pushing to optimum)
%   - Local Gradient:  ~ 200,000 (Resisting motion)
%   
%   Ratio ~ 10:1. The "Default Step (1.0)" faces a wall 10x steeper than 
%   the downward slope. It MUST get trapped.
%   The "Smart Step (5500)" sees the global drop of ~10^10 over its stride,
%   ignoring the local oscillation of ~10^5.
% -------------------------------------------------------------------------

% --- Problem Setup ---
n = 3;
x_opt = 100000 * ones(n, 1);
x0    = 110000 * ones(n, 1);

% Define Function: 
% Increased Amplitude 'A' from 4e4 to 2e5.
% Trap is now 5x deeper than before.
% OLD:
% fun = @(x) sum((x - x_opt).^2) + 2e5 * sum(cos(1.0 * x));

% NEW (Rigorous):
% Aligns the cosine valley exactly with the quadratic bottom.
fun = @(x) sum((x - x_opt).^2) + 2e5 * sum(1 - cos(x - x_opt));

fprintf('============================================================\n');
fprintf('TEST CASE: 3D Deep Trap Basin\n');
fprintf('Condition: Local Gradient (200k) >>> Global Gradient (20k)\n');
fprintf('============================================================\n');

%% 1. Run Default Strategy
opts_def = struct();
opts_def.num_blocks = n;
opts_def.alpha_init = 1.0; 
opts_def.iprint = 1;

fprintf('\n--- [1] Running Default Strategy (alpha = 1.0) ---\n');
[x_def, f_def, exit_def, out_def] = bds(fun, x0, opts_def);

%% 2. Run Smart Initialization Strategy

opts_new = struct();
opts_new.num_blocks = n;
opts_new.alpha_init = 'auto';
opts_new.iprint = 1;

fprintf('\n--- [2] Running Smart Init Strategy (Scaled Alpha) ---\n');
[x_new, f_new, exit_new, out_new] = bds(fun, x0, opts_new);

%% 3. Results Analysis

fprintf('\n============================================================\n');
fprintf('                  COMPARATIVE RESULTS (DEEP TRAP)           \n');
fprintf('============================================================\n');

dist_def = norm(x_def - x_opt);
dist_new = norm(x_new - x_opt);

% Threshold for "Global Basin" is looser due to high noise amplitude
status_def = "TRAPPED (Local Min)";
if dist_def < 5000, status_def = "SUCCESS (Global Basin)"; end

status_new = "TRAPPED (Local Min)";
if dist_new < 5000, status_new = "SUCCESS (Global Basin)"; end

fprintf('%-15s | %-10s | %-15s | %-s\n', 'Strategy', 'Func Evals', 'Final Dist', 'Result');
fprintf('%s\n', repmat('-', 1, 90));

fprintf('%-15s | %-10d | %-15.4e | %s\n', ...
    'Default (1.0)', out_def.funcCount, dist_def, status_def);

fprintf('%-15s | %-10d | %-15.4e | %s\n', ...
    'Smart Init', out_new.funcCount, dist_new, status_new);

fprintf('%s\n', repmat('-', 1, 90));

% Verification
if dist_def > 10000 && dist_new < 5000
    fprintf('\nCONCLUSION: SUCCESS. The deep trap caught the default strategy.\n');
else
    fprintf('\nCONCLUSION: Result is ambiguous. Trap might need frequency adjustment.\n');
end