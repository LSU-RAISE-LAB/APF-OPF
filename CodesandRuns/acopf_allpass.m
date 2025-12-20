function [results, success] = acopf_allpass(case_name, m, a_params, ipopt_overrides, dctype, loading_condition_per)

% ACOPF with rotated ALL-PASS rational kernel (YALMIP+IPOPT+DC pre-rotation)
% Usage:
%   [res, ok] = acopf_allpass('case118', 0);             % no RATE_A reduction
%   [res, ok] = acopf_allpass('case118', 10);            % 10% reduction on 20% of lines
%   [res, ok] = acopf_allpass('case118', 10, 0.52);      % K=1, specified a
%   [res, ok] = acopf_allpass('case118', 10, [0.3 0.2]); % K=2
%
% INPUTS
%   case_name        MATPOWER case (string) or mpc struct
%   m                Percentage reduction of RATE_A on 20% of in-service
%                    rated lines (those with smallest RATE_A).
%   a_params         (optional) vector of positive all-pass poles (K entries).
%                    If empty / omitted -> auto K=1 using band-based a = tan(Δ/2)/Δ.
%   ipopt_overrides  (optional) struct of fields to override IPOPT defaults:
%                    e.g., struct('max_iter',3000,'tol',1e-7)
%
% OUTPUTS
%   results          MATPOWER-like results struct:
%                    .success (1/0), .et (elapsed s), .f (objective),
%                    .bus, .gen, .branch (updated V, θ, Pg, Qg),
%                    .deltaE_rad, .deltaE_deg, .deltaE_from_bus, .deltaE_to_bus
%                    .cong_branch, .cong_bus
%   success          logical success flag (same as results.success)
%
% REQUIREMENTS
%   MATPOWER (loadcase, ext2int, int2ext, makeYbus, runpf/mpoption),
%   YALMIP, IPOPT.

if nargin < 1 || isempty(case_name), case_name = 'case118'; end
if nargin < 2 || isempty(m), m = 0; end
if nargin < 3, a_params = []; end
if nargin < 4, ipopt_overrides = struct(); end
if nargin < 5, dctype = 0; end
if nargin < 6 || isempty(loading_condition_per), loading_condition_per = 100; end
alpha = loading_condition_per/100;

yalmip('clear');
assert(exist('yalmip','file')==2,   'YALMIP not found on path.');
assert(exist('makeYbus','file')==2, 'MATPOWER not found on path.');
assert(exist('runpf','file')==2,    'MATPOWER runpf not found on path.');

% ---------- Load case & map to internal indexing ----------
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
% Bus
BUS_TYPE=2; PD=3; QD=4; VM=8; VA=9; VMAX=12; VMIN=13;
% Gen
GEN_BUS=1; PG=2; QG=3; QMAX=4; QMIN=5; PMAX=9; PMIN=10;
% Branch
F_BUS=1; T_BUS=2; BR_R=3; BR_X=4; BR_B=5; RATE_A=6; TAP=9; SHIFT=10; BR_STATUS=11; ANGMIN=12; ANGMAX=13;


% ---------- Load scaling ----------
bus(:,PD) = bus(:,PD) * alpha;
bus(:,QD) = bus(:,QD) * alpha;
mpc.bus   = bus;   % IMPORTANT: DC baseline uses mpc

% ---------- Basic data ----------
nb = size(bus,1); ng = size(gen,1);
Pd = bus(:,PD)/baseMVA; Qd = bus(:,QD)/baseMVA;
Vmin = bus(:,VMIN); Vmax = bus(:,VMAX);
ref_bus_idx = find(bus(:,BUS_TYPE)==3,1,'first'); assert(~isempty(ref_bus_idx),'No slack bus (type 3).');
genbus = gen(:,GEN_BUS);
Pgmin = gen(:,PMIN)/baseMVA; Pgmax = gen(:,PMAX)/baseMVA;
Qgmin = gen(:,QMIN)/baseMVA; Qgmax = gen(:,QMAX)/baseMVA;

% ---------- Network admittance & edge list for PF balance ----------
[Ybus,~,~] = makeYbus(mpc);           % complex
absY = abs(Ybus); phiY = angle(Ybus);
[bi, bj] = find(absY ~= 0);           % edges over nonzero Y_ij (internal indexing)
E = numel(bi);
absYe = zeros(E,1); phi_e = zeros(E,1);
for e = 1:E
    absYe(e) = absY(bi(e), bj(e));
    phi_e(e) = phiY(bi(e), bj(e));    % admittance angle φ_e
end
A = sparse(bi, (1:E)', 1, nb, E);     % rows sum outgoing edges at bus i

% ---------- Branch quantities for line flows ----------
on = (branch(:,BR_STATUS)==1);
on_idx  = find(on);                   % internal in-service branch indices

fbus_all  = branch(on,F_BUS);
tbus_all  = branch(on,T_BUS);
r         = branch(on,BR_R);
x         = branch(on,BR_X);
b_total   = branch(on,BR_B);
rateA     = branch(on,RATE_A);
tap_mag   = branch(on,TAP); tap_mag(tap_mag==0) = 1;
shift_deg = branch(on,SHIFT);
t         = tap_mag .* exp(1j*deg2rad(shift_deg));

% === NEW: apply m% reduction on 20% of in-service lines with smallest RATE_A ===
rateA_mod = rateA;
idx_candidates = find(~isnan(rateA_mod) & rateA_mod > 0);
if ~isempty(idx_candidates) && m ~= 0
    n_cand = numel(idx_candidates);
    n_select = max(1, round(0.90 * n_cand));  % 40% of rated, in-service lines
    [~, ord] = sort(rateA_mod(idx_candidates), 'ascend');  % smallest ratings first
    sel_local = ord(1:n_select);
    sel_idx = idx_candidates(sel_local);       % indices into "on" set
    rateA_mod(sel_idx) = rateA_mod(sel_idx) * (1 - m/100);
else
    sel_idx = [];
end
% === END NEW ===

y_series  = 1./(r + 1j*x);
y_sh      = 1j*(b_total);

a_f = (y_series + 0.5*y_sh)./(abs(t).^2);
b_f =  y_series ./ conj(t);
a_t = (y_series + 0.5*y_sh);
b_t =  y_series ./ t;

af_conj_r = real(conj(a_f)); af_conj_i = imag(conj(a_f));
at_conj_r = real(conj(a_t)); at_conj_i = imag(conj(a_t));
bf_mag    = abs(b_f);        ang_bf    = angle(b_f);
bt_mag    = abs(b_t);        ang_bt    = angle(b_t);

angmin = deg2rad(branch(on,ANGMIN));
angmax = deg2rad(branch(on,ANGMAX));
has_l  = ~isnan(angmin) & (angmin > -pi);
has_u  = ~isnan(angmax) & (angmax <  pi);

% ---------- Auto a_params from angle band if empty (equal poles, K=1) ----------
if isempty(a_params)
    if any(has_l | has_u)
        lo = angmin; lo(~has_l) = -pi/3;
        hi = angmax; hi(~has_u) =  pi/3;
        Delta = max(max(abs(lo)), max(abs(hi)));
    else
        Delta = pi/3; % fallback 60 deg
    end
    K = 1;
    a = tan(Delta/(2*K)) / max(Delta,1e-9);
    a_params = repmat(a,1,K);
end
a_params = a_params(:)';    % row
assert(all(a_params > 0), 'All-pass parameters must be > 0.');

% ---------- DC baseline for pre-rotation ----------
mpopt = mpoption('verbose',0,'out.all',0,'model','DC');
if dctype ==0 
res_dc = runpf(mpc, mpopt); tsolvedc = res_dc.et;
elseif dctype ==1
res_dc = rundcopf(mpc, mpopt); tsolvedc = res_dc.et;
end
assert(res_dc.success==1, 'DC PF failed; cannot build baseline.');
Va0_dc = deg2rad(res_dc.bus(:,VA));

% Baseline for PF-balance edges
deltaE0 = Va0_dc(bi) - Va0_dc(bj) - phi_e;
c0E = cos(deltaE0); s0E = sin(deltaE0);

% Baseline for line ends
d_f0 = Va0_dc(fbus_all) - Va0_dc(tbus_all) - ang_bf;
d_t0 = Va0_dc(tbus_all) - Va0_dc(fbus_all) - ang_bt;
cf0 = cos(d_f0); sf0 = sin(d_f0);
ct0 = cos(d_t0); st0 = sin(d_t0);

% ---------- Initial guesses ----------
Vm0 = max(min(bus(:,VM), Vmax), Vmin);
Va0 = bus(:,VA)*pi/180;
Pg0 = gen(:,PG)/baseMVA;
Qg0 = gen(:,QG)/baseMVA;

% ---------- IPOPT options ----------
ops = sdpsettings('solver','ipopt','usex0',1,'verbose',0,'cachesolvers',1, ...
                 'savesolveroutput',1,'savesolverinput',1);
try
    ops.ipopt.print_level = 0;
    ops.print_time = false;
    ops.ipopt.hessian_approximation = 'limited-memory';
    ops.ipopt.linear_solver = 'mumps';
    ops.ipopt.max_iter = 10000;
    ops.ipopt.tol = 1e-6;
    ops.ipopt.acceptable_tol = 1e-5;
    ops.ipopt.dual_inf_tol = 1e-6;
    ops.ipopt.constr_viol_tol = 1e-6;
catch
end
% Apply user overrides if provided
ipopt_fields = fieldnames(ipopt_overrides);
for ii=1:numel(ipopt_fields)
    try
        ops.ipopt.(ipopt_fields{ii}) = ipopt_overrides.(ipopt_fields{ii});
    catch
        % ignore invalid fields silently
    end
end

% ---------- Optimization model (ALL-PASS everywhere, rotated) ----------
V  = sdpvar(nb,1); th = sdpvar(nb,1);
Pg = sdpvar(ng,1); Qg = sdpvar(ng,1);

% PF balance via all-pass around DC baseline
resE = (th(bi) - th(bj) - phi_e) - deltaE0;      % residual vs baseline
[ReR, ImR] = allpass_cascade(resE, a_params);    % r(resE)
ReE = c0E .* ReR - s0E .* ImR;                   % approx cos(delta)
ImE = s0E .* ReR + c0E .* ImR;                   % approx sin(delta)
magE = V(bi).*V(bj).*absYe;
Pf_bus = A * (magE .* ReE);
Qf_bus = A * (magE .* ImE);
Pbal = Pf_bus; Qbal = Qf_bus;

% Line flows with rotated all-pass
rf = (th(fbus_all) - th(tbus_all) - ang_bf) - d_f0;
rt = (th(tbus_all) - th(fbus_all) - ang_bt) - d_t0;
[ReF, ImF] = allpass_cascade(rf, a_params);
[ReT, ImT] = allpass_cascade(rt, a_params);
cf = cf0 .* ReF - sf0 .* ImF;   sf = sf0 .* ReF + cf0 .* ImF;
ct = ct0 .* ReT - st0 .* ImT;   st = st0 .* ReT + ct0 .* ImT;

Vf = V(fbus_all); Vt = V(tbus_all);
Pf_br = af_conj_r .* (Vf.^2) - bf_mag .* (Vf.*Vt) .* cf;
Qf_br = af_conj_i .* (Vf.^2) - bf_mag .* (Vf.*Vt) .* sf;
Pt_br = at_conj_r .* (Vt.^2) - bt_mag .* (Vt.*Vf) .* ct;
Qt_br = at_conj_i .* (Vt.^2) - bt_mag .* (Vt.*Vf) .* st;

% Constraints
con = [];
con = [con, Vmin <= V <= Vmax, -pi <= th <= pi, th(ref_bus_idx)==0];
con = [con, Pgmin <= Pg <= Pgmax, Qgmin <= Qg <= Qgmax];

% Power balance using gen-bus incidence (handles multi-gen buses)
Ag = sparse(genbus, (1:ng)', 1, nb, ng);
con = [con, Ag*Pg - Pd == Pbal,  Ag*Qg - Qd == Qbal];

% Angle limits
if any(has_l), con = [con, th(fbus_all(has_l)) - th(tbus_all(has_l)) >= angmin(has_l)]; end
if any(has_u), con = [con, th(fbus_all(has_u)) - th(tbus_all(has_u)) <= angmax(has_u)]; end

% Thermal (apparent-power) limits where RATE_A present
Smax = rateA_mod/baseMVA;
idx_rate = find(~isnan(Smax) & Smax>0);
Smax_on = nan(size(rateA_mod));          % for diagnostics later
Smax_on(idx_rate) = Smax(idx_rate);

if ~isempty(idx_rate)
    con = [con, Pf_br(idx_rate).^2 + Qf_br(idx_rate).^2 <= Smax(idx_rate).^2];
    con = [con, Pt_br(idx_rate).^2 + Qt_br(idx_rate).^2 <= Smax(idx_rate).^2];
end

% Objective from gencost (type-2 only)
f = 0;
for kk = 1:ng
    gc = gencost(kk,:); if gc(1)~=2, continue; end
    ncost = gc(4); coeffs = gc(end-ncost+1:end);
    Pmw = Pg(kk)*baseMVA; tk = 0;
    for d = 1:ncost
        pow = ncost - d;
        tk = tk + coeffs(d) * Pmw^pow;
    end
    f = f + tk;
end

% Warm start from case only
assign(V, max(Vm0, Vmin + 1e-4));
assign(th, Va0); assign(Pg, Pg0); assign(Qg, Qg0);

% ---------- Solve ----------
sol = optimize(con, f, ops);
et = sol.solvertime;

ok = (sol.problem == 0);
if ~ok
    warning('acopf_allpass: %s', sol.info);
end



% ---------- Gather values ----------
Vm = value(V);
Va = value(th);                 % radians
Pg_pu = value(Pg); Qg_pu = value(Qg);
Pg_MW = Pg_pu * baseMVA; Qg_Mvar = Qg_pu * baseMVA;
obj = value(f);

% Compute mismatches for print diagnostics
Pg_bus = accumarray(genbus, Pg_pu, [nb,1], @sum, 0);
Qg_bus = accumarray(genbus, Qg_pu, [nb,1], @sum, 0);
Pmis = Pg_bus - Pd - value(Pbal);
Qmis = Qg_bus - Qd - value(Qbal);

% ---------- Line & bus congestion diagnostics (after solve) ----------
% Branch flows (numeric, physical trig)
Vf_val = Vm(fbus_all); Vt_val = Vm(tbus_all);
d_f_val = Va(fbus_all) - Va(tbus_all) - ang_bf;
d_t_val = Va(tbus_all) - Va(fbus_all) - ang_bt;
cf_val = cos(d_f_val); sf_val = sin(d_f_val);
ct_val = cos(d_t_val); st_val = sin(d_t_val);

Pf_val = af_conj_r.*(Vf_val.^2) - bf_mag.*(Vf_val.*Vt_val).*cf_val;
Qf_val = af_conj_i.*(Vf_val.^2) - bf_mag.*(Vf_val.*Vt_val).*sf_val;
Pt_val = at_conj_r.*(Vt_val.^2) - bt_mag.*(Vt_val.*Vf_val).*ct_val;
Qt_val = at_conj_i.*(Vt_val.^2) - bt_mag.*(Vt_val.*Vf_val).*st_val;

Sf_on = sqrt(Pf_val.^2 + Qf_val.^2);
St_on = sqrt(Pt_val.^2 + Qt_val.^2);

% Thermal congestion
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
% ---------- Angle difference per physical line ----------
dtheta_line_rad = ang_ft_val;               % radians
dtheta_line_deg = dtheta_line_rad * 180/pi; % degrees

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
% Using solved Va (radians) and phi_e.
% Direction convention: edge is from bus bj(e) -> bus bi(e),
% and deltaE = theta_to - theta_from - phi_e.
deltaE_all_rad = Va(bi) - Va(bj) - phi_e;

% Only keep off-diagonal Ybus entries (i ~= j)
is_offdiag = (bi ~= bj);
deltaE_rad = deltaE_all_rad(is_offdiag);
deltaE_deg = deltaE_rad * 180/pi;
bi_line   = bi(is_offdiag);   % "to" bus (internal)
bj_line   = bj(is_offdiag);   % "from" bus (internal)

% Map indices to EXTERNAL numbering when applicable
if had_order
    bus_map = mpc.order.bus.i2e;   % internal -> external
    deltaE_from_bus = bus_map(bj_line);
    deltaE_to_bus   = bus_map(bi_line);

    if isfield(mpc.order,'branch') && isfield(mpc.order.branch,'i2e')
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

% ---------- Populate MATPOWER-like results ----------
mpc_out = mpc; % internal indexing
mpc_out.bus(:,VM) = Vm;
mpc_out.bus(:,VA) = Va * 180/pi;
mpc_out.gen(:,PG) = Pg_MW;
mpc_out.gen(:,QG) = Qg_Mvar;

% Map back to external ordering if needed
if had_order
    results = int2ext(mpc_out);
else
    results = mpc_out;
end
% ---- Save full branch flows into results.branch (MATPOWER-style) ----
PF = 14; QF = 15; PT = 16; QT = 17;

% Initialize to NaN for all branches (including out-of-service)
results.branch(:, PF:QT) = NaN;

% Convert computed per-unit flows (for in-service "on" branches) to MW/MVAr
Pf_MW   = Pf_val * baseMVA;
Qf_MVAr = Qf_val * baseMVA;
Pt_MW   = Pt_val * baseMVA;
Qt_MVAr = Qt_val * baseMVA;

% Map internal in-service branch indices (on_idx) -> external branch indices
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


% Attach common fields
results.success = double(ok);
results.et      = et;
results.f       = obj;
results.Vm      = Vm;
results.Va      = Va * 180/pi;
results.Pg      = Pg_MW;
results.Qg      = Qg_Mvar;



% Try to extract iteration count from whatever the IPOPT interface returns
iters = NaN;
out = [];
if isfield(sol,'solveroutput'); out = sol.solveroutput; end

% Common patterns (depends on the IPOPT MATLAB interface you have)
if isstruct(out)
    if isfield(out,'iterations'), iters = out.iterations; end
    if isnan(iters) && isfield(out,'info') && isstruct(out.info) && isfield(out.info,'iter')
        iters = out.info.iter;
    end
    if isnan(iters) && isfield(out,'output') && isstruct(out.output) && isfield(out.output,'iterations')
        iters = out.output.iterations;
    end
end
results.ipopt_iter = iters;

% Attach deltaE diagnostics (all-pass version)
results.deltaE_rad      = deltaE_rad;
results.deltaE_deg      = deltaE_deg;
results.deltaE_from_bus = deltaE_from_bus;
results.deltaE_to_bus   = deltaE_to_bus;

% Attach congestion info
% Branch-level
results.cong_branch.branch_idx   = cong_branch_idx_ext(:);
results.cong_branch.from_bus     = results.branch(cong_branch_idx_ext, F_BUS);
results.cong_branch.to_bus       = results.branch(cong_branch_idx_ext, T_BUS);

Sf_cong   = Sf_on(cong_on_idx);
St_cong   = St_on(cong_on_idx);
Smax_cong = Smax_on(cong_on_idx);
if had_order
    bus_map = mpc.order.bus.i2e;

    dtheta_from_bus = bus_map(fbus_all);
    dtheta_to_bus   = bus_map(tbus_all);

    if isfield(mpc.order,'branch') && isfield(mpc.order.branch,'i2e')
        dtheta_branch_idx = mpc.order.branch.i2e(on_idx);
    else
        dtheta_branch_idx = on_idx;
    end
else
    dtheta_from_bus  = fbus_all;
    dtheta_to_bus    = tbus_all;
    dtheta_branch_idx = on_idx;
end
% ---------- Attach theta_i - theta_j diagnostics ----------
results.dtheta_line.branch_idx = dtheta_branch_idx(:);
results.dtheta_line.from_bus  = dtheta_from_bus(:);
results.dtheta_line.to_bus    = dtheta_to_bus(:);

results.dtheta_line.rad = dtheta_line_rad(:);
results.dtheta_line.deg = dtheta_line_deg(:);

results.cong_branch.Sf           = Sf_cong;
results.cong_branch.St           = St_cong;
results.cong_branch.Smax         = Smax_cong;
results.cong_branch.loading_from = Sf_cong ./ Smax_cong;
results.cong_branch.loading_to   = St_cong ./ Smax_cong;
results.cong_branch.thermal_flag = thermal_flag_on(cong_on_idx);
results.cong_branch.angle_flag   = angle_flag_on(cong_on_idx);
results.dctype           = dctype;
% Bus-level
results.cong_bus.bus_idx     = cong_bus_idx_ext(:);
results.cong_bus.Vm          = Vm(cong_bus_internal);
results.cong_bus.Va_deg      = Va(cong_bus_internal)*180/pi;
results.cong_bus.Pmis        = Pmis(cong_bus_internal);
results.cong_bus.Qmis        = Qmis(cong_bus_internal);
results.cong_bus.at_Vmin     = at_Vmin(cong_bus_internal);
results.cong_bus.at_Vmax     = at_Vmax(cong_bus_internal);
results.cong_bus.has_P_limit = Pbind_bus(cong_bus_internal);
results.cong_bus.has_Q_limit = Qbind_bus(cong_bus_internal);

% ---------- Print summary (MATPOWER-style spirit) ----------
if ischar(case_name) || (exist('isstring','builtin') && isstring(case_name))
    cname = char(case_name);
else
    cname = 'mpc';
end
fprintf('\n===== AC-OPF (ALL-PASS, DC pre-rotation) =====\n');
fprintf('Case: %s | K=%d | params: [%s]\n', cname, numel(a_params), num2str(a_params(:).',' %.4g'));
fprintf('Solver: YALMIP + IPOPT | hess=limited-memory | linear_solver=mumps\n');
fprintf('Status: %s\n', sol.info);
fprintf('Objective (cost): %.6f\n', obj);
fprintf('Elapsed time: %.3f s\n', et);
fprintf('IPOPT iterations: %d\n', iters);

fprintf('Pre-Rotation time: %.3f s\n', tsolvedc);
fprintf('Max |P-mis| (pu): %.3e | Max |Q-mis| (pu): %.3e\n', max(abs(Pmis)), max(abs(Qmis)));
fprintf('Vmag range (pu): [%.4f, %.4f]\n', min(Vm), max(Vm));
fprintf('Angle range (deg): [%.3f, %.3f]\n', min(Va)*180/pi, max(Va)*180/pi);
fprintf('Congested lines: %d\n', nCongLines);
fprintf('==============================================\n');

success = ok;
end

% ================== All-pass cascade helper ==================
function [ReR, ImR] = allpass_cascade(delta_vec, a_vec)
% r(δ) = Π_k (1 + j a_k δ) / (1 - j a_k δ), with |r(δ)|=1 for real δ
ReR = ones(size(delta_vec));
ImR = zeros(size(delta_vec));
for kk = 1:numel(a_vec)
    x  = a_vec(kk) .* delta_vec; % real vector
    nr = 1 - x.^2;               % Re{numerator}
    ni = 2 .* x;                 % Im{numerator}
    dr = 1 + x.^2;               % Re{denominator}
    newRe = (ReR .* nr - ImR .* ni) ./ dr;
    newIm = (ReR .* ni + ImR .* nr) ./ dr;
    ReR = newRe; ImR = newIm;
end
end
