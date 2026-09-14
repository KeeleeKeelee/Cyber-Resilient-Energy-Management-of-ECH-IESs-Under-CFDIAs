function model = dc_workload_flexibility_model()
%DC_WORKLOAD_FLEXIBILITY_MODEL  Spatiotemporal computing-workload flexibility
%   model of the data centres (DCs) used in the ECH-IES day-ahead dispatch.
%
%% ---------------------------------------------------------------- horizon
T  = 24;        % dispatch intervals
dt = 1;         % length of one dispatch interval

%% ------------------------------------------------------- DC nameplate data
%  Column 1 = DC1, column 2 = DC2.
N       = [800,  480];      % number of servers                    [-]
mu      = [150,  120];      % service rate of one server           [task/s]
U_max   = [0.90, 0.85];     % admissible server utilization        [-]
P_idle  = [150,  120];      % server power at zero utilization     [W]
P_peak  = [350,  280];      % server power at full utilization     [W]

n_dc = numel(N);

% Maximum computing rate of each DC
D_cap = N .* mu .* U_max;

%% ------------------------------------------------------- computing demands
% Interactive demand of the computing service, per interval.
Lambda_icw_total = 40000;                       % peak interactive demand [task/s]
icw_shape = [0.60 0.55 0.52 0.50 0.55 0.65 ...
             0.80 0.95 1.00 1.00 0.95 0.90 ...
             0.88 0.90 0.92 0.95 1.00 1.00 ...
             0.95 0.90 0.85 0.80 0.72 0.65];
Lambda_icw = Lambda_icw_total * icw_shape;

% Total number of batch tasks to be completed over the whole horizon.
Lambda_bcw = 5e9;

%% ------------------------------------------------ linear constraint blocks
% Variable ordering of the stacked workload vector z
nz    = 2 * n_dc * T;
blk   = 2 * n_dc;
i_icw = @(t) (t-1)*blk + (1:n_dc);
i_bcw = @(t) (t-1)*blk + n_dc + (1:n_dc);

% shared server-capacity envelope
A_cap = zeros(n_dc*T, nz);
b_cap = zeros(n_dc*T, 1);
for t = 1:T
    cols_icw = i_icw(t);
    cols_bcw = i_bcw(t);
    for i = 1:n_dc
        r = (t-1)*n_dc + i;
        A_cap(r, cols_icw(i)) = 1;
        A_cap(r, cols_bcw(i)) = 1;
        b_cap(r) = D_cap(i);
    end
end

% ICW conservation
Aeq_icw = zeros(T, nz);
beq_icw = zeros(T, 1);
for t = 1:T
    Aeq_icw(t, i_icw(t)) = 1;
    beq_icw(t) = Lambda_icw(t);
end

% BCW budget over the horizon
A_bcw = zeros(1, nz);
for t = 1:T
    A_bcw(i_bcw(t)) = -dt * 3600;
end
b_bcw = -Lambda_bcw;

%% --------------------------------------------------------------- assemble
model = struct();
model.T        = T;
model.dt       = dt;
model.n_dc     = n_dc;
model.N        = N;
model.mu       = mu;
model.U_max    = U_max;
model.P_idle   = P_idle;
model.P_peak   = P_peak;
model.D_cap    = D_cap;
model.Lambda_icw = Lambda_icw;
model.Lambda_bcw = Lambda_bcw;

model.Aineq    = [A_cap; A_bcw];
model.bineq    = [b_cap; b_bcw];
model.Aeq      = Aeq_icw;
model.beq      = beq_icw;
model.lb       = zeros(nz, 1);
model.index    = struct('icw', i_icw, 'bcw', i_bcw, 'nz', nz);

model.response = @response;

%% ------------------------------------------------------------- DC response
    function [P_IT, U] = response(D_icw, D_bcw)
        % Eqs. (3) and (2).  Workload in task/s, power returned in kW.
        U    = (D_icw + D_bcw) ./ (N(:) .* mu(:));
        P_IT = N(:) .* (P_idle(:) + (P_peak(:) - P_idle(:)) .* U) * 1e-3;
    end
end
