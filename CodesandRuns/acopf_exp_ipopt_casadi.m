function [results, success] = acopf_exp_ipopt_casadi(case_name, m, ipopt_overrides, loading_condition_per)

% AC-OPF with standard trigonometric kernel, NO YALMIP (CasADi + IPOPT).
% Usage:
%   [res, ok] = acopf_exp_ipopt_casadi('case118', 10);        % 10% reduction
%   [res, ok] = acopf_exp_ipopt_casadi('case118', 10, struct('tol',1e-7));
%
% Requires: MATPOWER on path; CasADi on path (import casadi.*).
%
% ADDED INPUT:
%   m                  - percentage reduction of RATE_A on 20% of in-service
%                        rated lines (those with smallest RATE_A).
%
% ADDED OUTPUTS (in results):
%   .deltaE_rad        - angle diffs δ_e = θ_i - θ_j - φ_e (radians) for off-diagonal Ybus entries
%   .deltaE_deg        - same in degrees
%   .deltaE_from_bus   - "from" bus index for each edge (external numbering)
%   .deltaE_to_bus     - "to" bus index for each edge (external numbering)
%
%   .cong_branch       - struct with congested line information:
%                        .branch_idx     (external branch row indices)
%                        .from_bus, .to_bus
%                        .Sf, .St        (|S| at from/to ends, pu)
%                        .Smax           (thermal rating, pu; NaN if none)
%                        .loading_from, .loading_to (Sf/Smax, St/Smax)
%                        .thermal_flag   (thermal limit binding)
%                        .angle_flag     (angle limit binding)
%
%   .cong_bus          - struct with congested bus information:
%                        .bus_idx        (external bus indices)
%                        .Vm, .Va_deg
%                        .Pmis, .Qmis    (pu)
%                        .at_Vmin, .at_Vmax
%                        .has_P_limit, .has_Q_limit

if nargin < 1 || isempty(case_name), case_name = 'case118'; end
if nargin < 2 || isempty(m), m = 0; end              % NEW: default no reduction
if nargin < 3, ipopt_overrides = struct(); end
if nargin < 4 || isempty(loading_condition_per), loading_condition_per = 100; end
alpha = loading_condition_per/100;

assert(exist('loadcase','file')==2, 'MATPOWER not on path.');
%assert(exist('casadi','file')==2,   'CasADi not on path.');

import casadi.*

% ---------- Load MATPOWER case & internal ordering ----------
mpc_in = loadcase(case_name);
had_order = isfield(mpc_in,'order');
if had_order
    if ~isfield(mpc_in.order,'state') || ~strcmp(mpc_in.order.state,'i')
        mpc = ext2int(mpc_in);
    else
        mpc = mpc_in;
    end
else
    mpc = ext2int(mpc_in);
end

baseMVA = mpc.baseMVA;
bus = mpc.bus; gen = mpc.gen; branch = mpc.branch; gencost = mpc.gencost;

% ---------- Indices ----------
BUS_TYPE=2; PD=3; QD=4; VM=8; VA=9; VMAX=12; VMIN=13;
GEN_BUS=1; PG=2; QG=3; QMAX=4; QMIN=5; PMAX=9; PMIN=10;
F_BUS=1; T_BUS=2; BR_R=3; BR_X=4; BR_B=5; RATE_A=6; TAP=9; SHIFT=10; BR_STATUS=11; ANGMIN=12; ANGMAX=13;
% ---------- Load scaling ----------
bus(:,PD) = bus(:,PD) * alpha;
bus(:,QD) = bus(:,QD) * alpha;
mpc.bus   = bus;

% ---------- Basic data ----------
nb = size(bus,1); ng = size(gen,1);
Pd = bus(:,PD)/baseMVA; Qd = bus(:,QD)/baseMVA;
Vmin = bus(:,VMIN); Vmax = bus(:,VMAX);
ref_bus_idx = find(bus(:,BUS_TYPE)==3,1,'first'); assert(~isempty(ref_bus_idx),'No slack bus.');
genbus = gen(:,GEN_BUS);
Pgmin = gen(:,PMIN)/baseMVA; Pgmax = gen(:,PMAX)/baseMVA;
Qgmin = gen(:,QMIN)/baseMVA; Qgmax = gen(:,QMAX)/baseMVA;

% ---------- Ybus for bus balances (include all nonzeros) ----------
[Ybus,~,~] = makeYbus(mpc);
absY = abs(Ybus); phiY = angle(Ybus);
[bi,bj] = find(absY~=0);                  % edges include diagonal (internal indexing)
E = numel(bi);
absYe = zeros(E,1); phi_e = zeros(E,1);
for e=1:E
    absYe(e)=absY(bi(e),bj(e));
    phi_e(e)=phiY(bi(e),bj(e));
end
A = sparse(bi,(1:E)',1,nb,E);             % accumulates per bus

% ---------- Branch data for thermal & angle limits ----------
on = (branch(:,BR_STATUS)==1);
on_idx = find(on);                         % internal branch indices for in-service lines

fbus_all = branch(on,F_BUS); tbus_all = branch(on,T_BUS);
r=branch(on,BR_R); x=branch(on,BR_X); b_total=branch(on,BR_B);
rateA=branch(on,RATE_A);
tap_mag=branch(on,TAP); tap_mag(tap_mag==0)=1;
shift_deg=branch(on,SHIFT);
t = tap_mag .* exp(1j*deg2rad(shift_deg));

% === NEW: apply m% reduction on 20% of in-service lines with smallest RATE_A ===
rateA_mod = rateA;
idx_candidates = find(~isnan(rateA_mod) & rateA_mod>0);
if ~isempty(idx_candidates) && m~=0
    n_cand = numel(idx_candidates);
    n_select = max(1, round(0.90 * n_cand));   % 40% of the rated, in-service lines
    [~, ord] = sort(rateA_mod(idx_candidates), 'ascend');  % smallest ratings = more likely congested
    sel_local = ord(1:n_select);
    sel_idx = idx_candidates(sel_local);       % indices into on-branch set
    rateA_mod(sel_idx) = rateA_mod(sel_idx) * (1 - m/100);
else
    sel_idx = [];
end
% === END NEW ===

y_series = 1./(r+1j*x);
y_sh     = 1j*(b_total);

a_f=(y_series+0.5*y_sh)./(abs(t).^2);
b_f= y_series./conj(t);
a_t=(y_series+0.5*y_sh);
b_t= y_series./t;

af_conj_r=real(conj(a_f)); af_conj_i=imag(conj(a_f));
at_conj_r=real(conj(a_t)); at_conj_i=imag(conj(a_t));
bf_mag=abs(b_f);           ang_bf=angle(b_f);
bt_mag=abs(b_t);           ang_bt=angle(b_t);

angmin=deg2rad(branch(on,ANGMIN)); angmax=deg2rad(branch(on,ANGMAX));
has_l=~isnan(angmin)&(angmin>-pi);
has_u=~isnan(angmax)&(angmax< pi);

% ---------- Variables (CasADi SX) ----------
V  = SX.sym('V',  nb);      % voltage magnitudes
th = SX.sym('th', nb);      % angles (rad)
Pg = SX.sym('Pg', ng);      % active power (pu)
Qg = SX.sym('Qg', ng);      % reactive power (pu)

% ---------- Initial guesses ----------
Vm0 = max(min(bus(:,VM),Vmax),Vmin);
Va0 = bus(:,VA)*pi/180;
Pg0 = gen(:,PG)/baseMVA;
Qg0 = gen(:,QG)/baseMVA;

x  = [V; th; Pg; Qg];
x0 = [Vm0; Va0; Pg0; Qg0];
lbx = [Vmin; -pi*ones(nb,1); Pgmin; Qgmin];
ubx = [Vmax;  pi*ones(nb,1); Pgmax; Qgmax];
% Fix slack angle = 0
lbx(nb + ref_bus_idx) = 0;
ubx(nb + ref_bus_idx) = 0;

% ---------- Bus power balance (trig kernel) ----------
deltaE = th(bi) - th(bj) - phi_e;
ReE = cos(deltaE);
ImE = sin(deltaE);
magE = V(bi).*V(bj).*absYe;

Pf_bus = A * (magE .* ReE);
Qf_bus = A * (magE .* ImE);

Ag = sparse(genbus,(1:ng)',1,nb,ng);
Pg_bus = Ag * Pg;
Qg_bus = Ag * Qg;

gP = Pg_bus - Pd - Pf_bus;            % = 0
gQ = Qg_bus - Qd - Qf_bus;            % = 0

% ---------- Angle limits ----------
ineq = SX([]); lb_ineq = []; ub_ineq = [];
if any(has_l)
    expr = th(fbus_all(has_l)) - th(tbus_all(has_l)) - angmin(has_l);
    ineq = [ineq; expr];
    lb_ineq = [lb_ineq; zeros(nnz(has_l),1)];
    ub_ineq = [ub_ineq;  inf(nnz(has_l),1)];
end
if any(has_u)
    expr = angmax(has_u) - (th(fbus_all(has_u)) - th(tbus_all(has_u)));
    ineq = [ineq; expr];
    lb_ineq = [lb_ineq; zeros(nnz(has_u),1)];
    ub_ineq = [ub_ineq;  inf(nnz(has_u),1)];
end

% ---------- Thermal (apparent power) limits ----------
Vf = V(fbus_all); Vt = V(tbus_all);
d_f = th(fbus_all)-th(tbus_all)-ang_bf;
d_t = th(tbus_all)-th(fbus_all)-ang_bt;
cf = cos(d_f); sf = sin(d_f);
ct = cos(d_t); st = sin(d_t);

Pf_br = af_conj_r.*(Vf.^2) - bf_mag.*(Vf.*Vt).*cf;
Qf_br = af_conj_i.*(Vf.^2) - bf_mag.*(Vf.*Vt).*sf;
Pt_br = at_conj_r.*(Vt.^2) - bt_mag.*(Vt.*Vf).*ct;
Qt_br = at_conj_i.*(Vt.^2) - bt_mag.*(Vt.*Vf).*st;

% NEW: use modified rateA_mod for Smax
Smax = rateA_mod/baseMVA;
idx_rate = find(~isnan(Smax) & Smax>0);
Smax_on = nan(size(rateA_mod));              % for later diagnostics
Smax_on(idx_rate) = Smax(idx_rate);

if ~isempty(idx_rate)
    Sf2 = Pf_br(idx_rate).^2 + Qf_br(idx_rate).^2;
    St2 = Pt_br(idx_rate).^2 + Qt_br(idx_rate).^2;
    % Inequalities: Smax^2 - S^2 >= 0
    ineq = [ineq;  (Smax(idx_rate).^2 - Sf2); (Smax(idx_rate).^2 - St2)];
    lb_ineq = [lb_ineq; zeros(2*numel(idx_rate),1)];
    ub_ineq = [ub_ineq;  inf(2*numel(idx_rate),1)];
end

% ---------- Stack constraints ----------
g   = [gP; gQ; ineq];
lbg = [zeros(2*nb,1); lb_ineq];
ubg = [zeros(2*nb,1); ub_ineq];

% ---------- Objective (MATPOWER type-2 poly) ----------
f = SX(0);
for kk=1:ng
    gc = gencost(kk,:); if gc(1)~=2, continue; end
    ncost = gc(4);
    coeffs = gc(end-ncost+1:end);
    Pmw = Pg(kk)*baseMVA;
    tk = 0;
    for d=1:ncost
        pow = ncost - d;
        tk = tk + coeffs(d)*Pmw^pow;
    end
    f = f + tk;
end

% ---------- NLP + IPOPT options ----------
nlp = struct('x',x,'f',f,'g',g);
opts = struct;
opts.ipopt.print_level = 0;
opts.print_time = false;
opts.ipopt.hessian_approximation = 'limited-memory';
opts.ipopt.linear_solver = 'mumps';
opts.ipopt.max_iter = 10000;
opts.ipopt.tol = 1e-6;
opts.ipopt.acceptable_tol = 1e-5;
opts.ipopt.dual_inf_tol = 1e-6;
opts.ipopt.constr_viol_tol = 1e-6;



% user overrides
fn = fieldnames(ipopt_overrides);
for k=1:numel(fn)
    opts.ipopt.(fn{k}) = ipopt_overrides.(fn{k});
end

solver = nlpsol('solver','ipopt',nlp,opts);

% ---------- Solve ----------
t0 = tic;
sol = solver('x0',x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg);
et = toc(t0);

stats = solver.stats;    % CasADi statistics struct

% Depending on CasADi version, stats may be nested (nlpsol.*) or flat.
if isfield(stats, 'nlpsol')
    iter_count = stats.nlpsol.iter_count;
    iters      = stats.nlpsol.iterations;      % struct with fields obj, inf_pr, ...
else
    iter_count = stats.iter_count;
    iters      = stats.iterations;
end

obj_hist = iters.obj(:);   % column vector of objective values per iteration


xsol = full(sol.x);
Vm = xsol(1:nb);
Va = xsol(nb+(1:nb));
Pg_pu = xsol(2*nb+(1:ng));
Qg_pu = xsol(2*nb+ng+(1:ng));

Pg_MW = Pg_pu*baseMVA; Qg_Mvar = Qg_pu*baseMVA;
obj = full(sol.f);

% ---------- Diagnostics (mismatch) ----------
deltaE_val = Va(bi) - Va(bj) - phi_e;
Pf_bus_val = A * (Vm(bi).*Vm(bj).*absYe .* cos(deltaE_val));
Qf_bus_val = A * (Vm(bi).*Vm(bj).*absYe .* sin(deltaE_val));
Pg_bus_val = accumarray(genbus, Pg_pu, [nb,1], @sum, 0);
Qg_bus_val = accumarray(genbus, Qg_pu, [nb,1], @sum, 0);
Pmis = Pg_bus_val - Pd - Pf_bus_val;
Qmis = Qg_bus_val - Qd - Qf_bus_val;

% ---------- Line & bus congestion diagnostics (after solve) ----------
% Branch flows (numeric)
Vf_val = Vm(fbus_all); Vt_val = Vm(tbus_all);
d_f_val = Va(fbus_all)-Va(tbus_all)-ang_bf;
d_t_val = Va(tbus_all)-Va(fbus_all)-ang_bt;
cf_val = cos(d_f_val); sf_val = sin(d_f_val);
ct_val = cos(d_t_val); st_val = sin(d_t_val);

Pf_val = af_conj_r.*(Vf_val.^2) - bf_mag.*(Vf_val.*Vt_val).*cf_val;
Qf_val = af_conj_i.*(Vf_val.^2) - bf_mag.*(Vf_val.*Vt_val).*sf_val;
Pt_val = at_conj_r.*(Vt_val.^2) - bt_mag.*(Vt_val.*Vf_val).*ct_val;
Qt_val = at_conj_i.*(Vt_val.^2) - bt_mag.*(Vt_val.*Vf_val).*st_val;

Sf_on = sqrt(Pf_val.^2 + Qf_val.^2);
St_on = sqrt(Pt_val.^2 + Qt_val.^2);

% Thermal congestion (relative to Smax_on)
tol_loading = 1e-3;
thermal_flag_on = false(size(fbus_all));
if ~isempty(idx_rate)
    Sf_idx = Sf_on(idx_rate);
    St_idx = St_on(idx_rate);
    ratio_f = Sf_idx ./ Smax(idx_rate);
    ratio_t = St_idx ./ Smax(idx_rate);
    thermal_cong_idx = (ratio_f >= (1 - tol_loading)) | (ratio_t >= (1 - tol_loading));
    thermal_flag_on(idx_rate(thermal_cong_idx)) = true;
end

% Angle congestion
tol_ang = 1e-4;  % rad
ang_ft_val = Va(fbus_all) - Va(tbus_all);
angle_cong_l = false(size(fbus_all));
angle_cong_u = false(size(fbus_all));
if any(has_l)
    angle_cong_l(has_l) = abs(ang_ft_val(has_l) - angmin(has_l)) <= tol_ang;
end
if any(has_u)
    angle_cong_u(has_u) = abs(ang_ft_val(has_u) - angmax(has_u)) <= tol_ang;
end
angle_flag_on = angle_cong_l | angle_cong_u;

% Lines congested by either thermal or angle limits
line_flag_on = thermal_flag_on | angle_flag_on;
cong_on_idx = find(line_flag_on);
nCongLines = numel(cong_on_idx);

% Map congested on-lines to global internal branch indices
cong_branch_internal = on_idx(cong_on_idx);

% ---------- Angle difference per physical line (AFTER solve) ----------
theta_f = Va(fbus_all);      % rad
theta_t = Va(tbus_all);      % rad

dtheta_line_rad = theta_f - theta_t;
dtheta_line_deg = dtheta_line_rad * 180/pi;

if had_order
    bus_map = mpc.order.bus.i2e;
    dtheta_from_bus = bus_map(fbus_all);
    dtheta_to_bus   = bus_map(tbus_all);

    if isfield(mpc.order, 'branch') && isfield(mpc.order.branch,'i2e')
        dtheta_branch_idx = mpc.order.branch.i2e(on_idx);
    else
        dtheta_branch_idx = on_idx;
    end
else
    dtheta_from_bus = fbus_all;
    dtheta_to_bus   = tbus_all;
    dtheta_branch_idx = on_idx;
end


% Bus congestion:
tol_V = 1e-4;
tol_PQ = 1e-5;
at_Vmin = abs(Vm - Vmin) <= tol_V;
at_Vmax = abs(Vm - Vmax) <= tol_V;

at_Pmin_gen = abs(Pg_pu - Pgmin) <= tol_PQ;
at_Pmax_gen = abs(Pg_pu - Pgmax) <= tol_PQ;
at_Qmin_gen = abs(Qg_pu - Qgmin) <= tol_PQ;
at_Qmax_gen = abs(Qg_pu - Qgmax) <= tol_PQ;

Pbind_bus = accumarray(genbus, double(at_Pmin_gen | at_Pmax_gen), [nb,1], @max, 0) > 0;
Qbind_bus = accumarray(genbus, double(at_Qmin_gen | at_Qmax_gen), [nb,1], @max, 0) > 0;

bus_flag = at_Vmin | at_Vmax | Pbind_bus | Qbind_bus;
cong_bus_internal = find(bus_flag);

% ---------- Compute deltaE outputs (AFTER solve) ----------
% Direction convention: edge is from bus bj(e) -> bus bi(e),
% and deltaE = theta_to - theta_from - phi_e.
is_offdiag = (bi ~= bj);
deltaE_rad = deltaE_val(is_offdiag);
deltaE_deg = deltaE_rad * 180/pi;
bi_line    = bi(is_offdiag);   % "to" bus (internal)
bj_line    = bj(is_offdiag);   % "from" bus (internal)

% Map indices to EXTERNAL numbering when applicable
if had_order
    bus_map = mpc.order.bus.i2e;   % internal -> external
    deltaE_from_bus = bus_map(bj_line);
    deltaE_to_bus   = bus_map(bi_line);

    if isfield(mpc.order, 'branch') && isfield(mpc.order.branch, 'i2e')
        br_map = mpc.order.branch.i2e;
        cong_branch_idx_ext = br_map(cong_branch_internal);
    else
        cong_branch_idx_ext = cong_branch_internal;
    end

    cong_bus_idx_ext = bus_map(cong_bus_internal);
else
    deltaE_from_bus = bj_line;
    deltaE_to_bus   = bi_line;
    cong_branch_idx_ext = cong_branch_internal;
    cong_bus_idx_ext = cong_bus_internal;
end

% ---------- Pack results ----------
mpc_out = mpc;
mpc_out.bus(:,VM) = Vm;
mpc_out.bus(:,VA) = Va*180/pi;
mpc_out.gen(:,PG) = Pg_MW;
mpc_out.gen(:,QG) = Qg_Mvar;

if had_order
    results = int2ext(mpc_out);
else
    results = mpc_out;
end
rs = solver.stats.return_status;
ok = strcmp(rs,'Solve_Succeeded') || strcmp(rs,'Solved_To_Acceptable_Level');
results.ipopt_iter      = iter_count;
results.ipopt_obj_hist  = obj_hist;

% ---- Save full branch flows into results.branch (MATPOWER-style) ----
PF = 14; QF = 15; PT = 16; QT = 17;

% Ensure columns exist; initialize to NaN for all branches (incl out-of-service)
results.branch(:, PF:QT) = NaN;

% Convert per-unit flows (computed for in-service "on" branches) to MW/MVAr
Pf_MW   = Pf_val * baseMVA;
Qf_MVAr = Qf_val * baseMVA;
Pt_MW   = Pt_val * baseMVA;
Qt_MVAr = Qt_val * baseMVA;

% Map internal in-service branch indices -> external branch indices
if had_order
    if isfield(mpc.order,'branch') && isfield(mpc.order.branch,'i2e')
        br_map = mpc.order.branch.i2e;     % internal -> external
        ext_on_idx = br_map(on_idx);       % external indices for in-service branches
    else
        ext_on_idx = on_idx;
    end
else
    ext_on_idx = on_idx; % internal == external
end

% Write flows for ALL in-service branches
results.branch(ext_on_idx, PF) = Pf_MW;
results.branch(ext_on_idx, QF) = Qf_MVAr;
results.branch(ext_on_idx, PT) = Pt_MW;
results.branch(ext_on_idx, QT) = Qt_MVAr;



results.success = double(ok);
results.et      = et;
results.f       = obj;
results.Vm = Vm;
results.Va = Va * 180/pi;
results.Pg = Pg_MW;
results.Qg = Qg_Mvar;

% Attach deltaE diagnostics
results.deltaE_rad      = deltaE_rad;
results.deltaE_deg      = deltaE_deg;
results.deltaE_from_bus = deltaE_from_bus;
results.deltaE_to_bus   = deltaE_to_bus;

% Attach congestion info
% Branch-level
results.cong_branch.branch_idx   = cong_branch_idx_ext(:);
results.cong_branch.from_bus     = results.branch(cong_branch_idx_ext, F_BUS);
results.cong_branch.to_bus       = results.branch(cong_branch_idx_ext, T_BUS);

% create Smax_on in external order by mapping through "on" mask
Sf_cong = Sf_on(cong_on_idx);
St_cong = St_on(cong_on_idx);
Smax_cong = Smax_on(cong_on_idx);

results.cong_branch.Sf          = Sf_cong;
results.cong_branch.St          = St_cong;
results.cong_branch.Smax        = Smax_cong;
results.cong_branch.loading_from = Sf_cong ./ Smax_cong;
results.cong_branch.loading_to   = St_cong ./ Smax_cong;
results.cong_branch.thermal_flag = thermal_flag_on(cong_on_idx);
results.cong_branch.angle_flag   = angle_flag_on(cong_on_idx);

% Bus-level
results.cong_bus.bus_idx    = cong_bus_idx_ext(:);
results.cong_bus.Vm         = Vm(cong_bus_internal);
results.cong_bus.Va_deg     = Va(cong_bus_internal)*180/pi;
results.cong_bus.Pmis       = Pmis(cong_bus_internal);
results.cong_bus.Qmis       = Qmis(cong_bus_internal);
results.cong_bus.at_Vmin    = at_Vmin(cong_bus_internal);
results.cong_bus.at_Vmax    = at_Vmax(cong_bus_internal);
results.cong_bus.has_P_limit = Pbind_bus(cong_bus_internal);
results.cong_bus.has_Q_limit = Qbind_bus(cong_bus_internal);

results.dtheta_line.branch_idx = dtheta_branch_idx(:);
results.dtheta_line.from_bus  = dtheta_from_bus(:);
results.dtheta_line.to_bus    = dtheta_to_bus(:);

results.dtheta_line.rad = dtheta_line_rad(:);
results.dtheta_line.deg = dtheta_line_deg(:);


% ---------- Print summary ----------
if ischar(case_name) || (exist('isstring','builtin') && isstring(case_name))
    cname = char(case_name);
else
    cname = 'mpc';
end
fprintf('\n===== AC-OPF (EXP baseline, CasADi+IPOPT) =====\n');
fprintf('Case: %s\n', cname);
fprintf('Solver: CasADi + IPOPT (exp kernel)\n');
fprintf('Status: %s\n', rs);
fprintf('Objective (cost): %.6f\n', obj);
fprintf('Elapsed time: %.3f s\n', et);
fprintf('IPOPT iterations: %d\n', iter_count);

fprintf('Max |P-mis| (pu): %.3e | Max |Q-mis| (pu): %.3e\n', max(abs(Pmis)), max(abs(Qmis)));
fprintf('Vmag range (pu): [%.4f, %.4f]\n', min(Vm), max(Vm));
fprintf('Angle range (deg): [%.3f, %.3f]\n', min(Va)*180/pi, max(Va)*180/pi);
fprintf('Congested lines: %d\n', nCongLines);   % NEW summary line
fprintf('===============================================\n');

success = ok;


try
    solver.delete();
catch
end
clear solver nlp x g f

end
