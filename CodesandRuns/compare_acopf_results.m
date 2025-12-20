function report = compare_acopf_results(res1, ok1, res2, ok2, varargin)
%COMPARE_ACOPF_RESULTS Clean comparison: EXP vs ALL-PASS
% Outputs:
%   report.results.exp = res1
%   report.results.ap  = res2
%   report.errors      = max/min/avg for Pg,Qg,Vm,Va (absolute)
%   report.objective   = f1,f2,df,gap_percent
%   report.stats       = dtheta_line & deltaE stats (from res structs)
%   report.congestion  = branch-only congestion (from res structs)
%   report.feasibility = res2 checked under TRUE trig AC constraints (bus balance, V, gen, thermal, angles)
%
% NEW:
%   Separate tolerances for feasibility checks:
%     tolP   - bus P balance (pu)
%     tolQ   - bus Q balance (pu)
%     tolV   - voltage magnitude (pu)
%     tolPg  - generator P limits (pu)
%     tolQg  - generator Q limits (pu)
%     tolAng - angle-difference limits (rad)
%     tolTh  - thermal limit exceedance (pu)
%
% Backward compatible:
%   'tol' still accepted and acts as the default for all tolerances,
%   unless a specific tolerance (tolP, tolV, ...) is provided.

% ---------------- Parse options ----------------
p = inputParser;

% Backward-compatible global tol (default)
p.addParameter('tol',      1e-1, @(x) isnumeric(x) && isscalar(x) && x>0);

% NEW: per-constraint tolerances (optional; if empty -> inherit from 'tol')
p.addParameter('tolP',   [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
p.addParameter('tolQ',   [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
p.addParameter('tolV',   [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
p.addParameter('tolPg',  [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
p.addParameter('tolQg',  [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
p.addParameter('tolAng', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
p.addParameter('tolTh',  [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));

p.addParameter('list_top', 10,   @(x) isnumeric(x) && isscalar(x) && x>=0);
p.parse(varargin{:});

tol_global = p.Results.tol;

% inherit any empty per-constraint tol from global tol
tolP   = inherit_tol(p.Results.tolP,   tol_global);
tolQ   = inherit_tol(p.Results.tolQ,   tol_global);
tolV   = inherit_tol(p.Results.tolV,   tol_global);
tolPg  = inherit_tol(p.Results.tolPg,  tol_global);
tolQg  = inherit_tol(p.Results.tolQg,  tol_global);
tolAng = inherit_tol(p.Results.tolAng, tol_global);
tolTh  = inherit_tol(p.Results.tolTh,  tol_global);

list_top = p.Results.list_top;

% ---------------- Indices ----------------
BUS_TYPE=2; PD=3; QD=4; VM=8; VA=9; VMAX=12; VMIN=13;
GEN_BUS=1; PG=2; QG=3; QMAX=4; QMIN=5; PMAX=9; PMIN=10;
F_BUS=1; T_BUS=2; BR_R=3; BR_X=4; BR_B=5; RATE_A=6; TAP=9; SHIFT=10; BR_STATUS=11; ANGMIN=12; ANGMAX=13;

% ---------------- Output struct (intact results) ----------------
report = struct();
report.results = struct();
report.results.exp = res1;
report.results.ap  = res2;

% ---------------- Variable errors (absolute) ----------------
vfields = {'Pg','Qg','Vm','Va'};
for k=1:numel(vfields)
    fn = vfields{k};
    if ~isfield(res1, fn) || ~isfield(res2, fn) || isempty(res1.(fn)) || isempty(res2.(fn))
        report.errors.(fn) = struct('max',NaN,'min',NaN,'avg',NaN);
    else
        d = abs(res2.(fn)(:) - res1.(fn)(:));
        report.errors.(fn) = struct('max',max(d),'min',min(d),'avg',mean(d));
    end
end

% ---------------- Branch flow errors (absolute) ----------------
PF = 14; QF = 15; PT = 16; QT = 17;

flow_fields = {'PF','QF','PT','QT'};
flow_cols   = [ PF , QF , PT , QT ];

for k = 1:numel(flow_fields)
    fn = flow_fields{k};
    c  = flow_cols(k);

    if ~isfield(res1,'branch') || ~isfield(res2,'branch') || ...
       size(res1.branch,2) < c || size(res2.branch,2) < c
        report.errors.(fn) = struct('max',NaN,'min',NaN,'avg',NaN);
        continue;
    end

    x1 = res1.branch(:,c);
    x2 = res2.branch(:,c);

    % Use only entries where BOTH are finite (skips out-of-service NaNs)
    mask = isfinite(x1) & isfinite(x2);

    if ~any(mask)
        report.errors.(fn) = struct('max',NaN,'min',NaN,'avg',NaN);
    else
        d = abs(x2(mask) - x1(mask));
        report.errors.(fn) = struct('max',max(d),'min',min(d),'avg',mean(d));
    end
end


% ---------------- Objective gap ----------------
assert(isfield(res1,'f') && isfield(res2,'f'), 'res1/res2 must contain field .f');
f1 = res1.f; f2 = res2.f;
report.objective.f_exp = f1;
report.objective.f_ap  = f2;
report.objective.df    = (f2 - f1);
report.objective.gap_percent = 100 * abs(f2 - f1) / max(1e-12, abs(f1));

% ---------------- Timing & iterations ----------------
report.timing = struct();
report.timing.exp.et_s = get_scalar_or_nan(res1, 'et');
report.timing.ap.et_s  = get_scalar_or_nan(res2, 'et');
report.timing.exp.ipopt_iter = get_scalar_or_nan(res1, 'ipopt_iter');
report.timing.ap.ipopt_iter  = get_scalar_or_nan(res2, 'ipopt_iter');

te = report.timing.exp.et_s;
ta = report.timing.ap.et_s;
if ~isnan(te) && ~isnan(ta) && te > 0
    report.timing.delta_s = ta - te;
    report.timing.delta_percent = 100*(ta - te)/te;
    report.timing.speedup = te / ta;
else
    report.timing.delta_s = NaN;
    report.timing.delta_percent = NaN;
    report.timing.speedup = NaN;
end

% ---------------- Stats from solver outputs (direct) ----------------
report.stats = struct();
report.stats.dtheta = stats_from_vector( ...
    getfield_safe(res1, {'dtheta_line','deg'}), ...
    getfield_safe(res2, {'dtheta_line','deg'}) );
report.stats.deltaE = stats_from_vector( ...
    getfield_safe(res1, {'deltaE_deg'}), ...
    getfield_safe(res2, {'deltaE_deg'}) );

% ---------------- Branch congestion (branch-only) ----------------
report.congestion = struct();
report.congestion.exp = extract_branch_congestion(res1);
report.congestion.ap  = extract_branch_congestion(res2);
report.congestion.match = compare_branch_congestion(report.congestion.exp, report.congestion.ap);

% ---------------- TRUE TRIG feasibility check using res2 ----------------
mpc2 = to_internal_copy(res2);
baseMVA = mpc2.baseMVA;
bus = mpc2.bus; gen = mpc2.gen; br = mpc2.branch;

nb = size(bus,1);
ng = size(gen,1);

Pd = bus(:,PD)/baseMVA;
Qd = bus(:,QD)/baseMVA;

Vmin = bus(:,VMIN); Vmax = bus(:,VMAX);
Pgmin = gen(:,PMIN)/baseMVA; Pgmax = gen(:,PMAX)/baseMVA;
Qgmin = gen(:,QMIN)/baseMVA; Qgmax = gen(:,QMAX)/baseMVA;

on = (br(:,BR_STATUS)==1);
fbus = br(on,F_BUS); tbus = br(on,T_BUS);

% Build Ybus edge list (includes diagonal)
[Ybus,~,~] = makeYbus(mpc2);
absY = abs(Ybus); phiY = angle(Ybus);
[bi,bj] = find(absY ~= 0);
E = numel(bi);
absYe = zeros(E,1); phi_e = zeros(E,1);
for e=1:E
    absYe(e)=absY(bi(e),bj(e));
    phi_e(e)=phiY(bi(e),bj(e));
end
A = sparse(bi,(1:E)',1,nb,E);

genbus = gen(:,GEN_BUS);
Ag = sparse(genbus,(1:ng)',1,nb,ng);

% Solution values from res2 (res2.Va is degrees)
Vm = res2.Vm(:);
Va = deg2rad(res2.Va(:));
Pg_pu_gen = res2.Pg(:)/baseMVA;
Qg_pu_gen = res2.Qg(:)/baseMVA;

Pg_bus = accumarray(genbus, Pg_pu_gen, [nb,1], @sum, 0);
Qg_bus = accumarray(genbus, Qg_pu_gen, [nb,1], @sum, 0);


% ================== EXP vs AP elementwise "error violations" (pu/rad) ==================
% Uses same tolerance set as feasibility:
%   tolP, tolQ  : branch flows (PF/PT, QF/QT) in pu
%   tolV        : Vm in pu
%   tolPg, tolQg: generator Pg/Qg in pu
%   tolAng      : angles in rad (Va diff)

report.errcheck = struct();
report.errcheck.tols = struct('tolP_pu',tolP,'tolQ_pu',tolQ,'tolV_pu',tolV,'tolPg_pu',tolPg,'tolQg_pu',tolQg,'tolAng_rad',tolAng);

% --- Generator injections: convert MW/MVAr -> pu
dPg_pu = abs(res2.Pg(:) - res1.Pg(:)) / baseMVA;
dQg_pu = abs(res2.Qg(:) - res1.Qg(:)) / baseMVA;

viol_Pg = dPg_pu > tolPg;
viol_Qg = dQg_pu > tolQg;

% --- Voltage magnitude already in pu
dVm = abs(res2.Vm(:) - res1.Vm(:));
viol_Vm = dVm > tolV;

% --- Bus angles: degrees -> rad
dVa_rad = abs(deg2rad(res2.Va(:)) - deg2rad(res1.Va(:)));
viol_Va = dVa_rad > tolAng;

% --- Branch flows: compare PF/PT together for P, QF/QT together for Q (pu)
PF = 14; QF = 15; PT = 16; QT = 17;

% Safety: allow NaNs/out-of-service -> ignore those entries
x1PF = res1.branch(:,PF); x2PF = res2.branch(:,PF);
x1PT = res1.branch(:,PT); x2PT = res2.branch(:,PT);
x1QF = res1.branch(:,QF); x2QF = res2.branch(:,QF);
x1QT = res1.branch(:,QT); x2QT = res2.branch(:,QT);

mPF = isfinite(x1PF) & isfinite(x2PF);
mPT = isfinite(x1PT) & isfinite(x2PT);
mQF = isfinite(x1QF) & isfinite(x2QF);
mQT = isfinite(x1QT) & isfinite(x2QT);

dPF_pu = NaN(size(x1PF)); dPT_pu = NaN(size(x1PT));
dQF_pu = NaN(size(x1QF)); dQT_pu = NaN(size(x1QT));

dPF_pu(mPF) = abs(x2PF(mPF) - x1PF(mPF)) / baseMVA;
dPT_pu(mPT) = abs(x2PT(mPT) - x1PT(mPT)) / baseMVA;
dQF_pu(mQF) = abs(x2QF(mQF) - x1QF(mQF)) / baseMVA;
dQT_pu(mQT) = abs(x2QT(mQT) - x1QT(mQT)) / baseMVA;

% Combine PF/PT into one "P-flow mismatch" flag per physical branch row
viol_Pflow = false(size(x1PF));
viol_Qflow = false(size(x1QF));

viol_Pflow(mPF) = viol_Pflow(mPF) | (dPF_pu(mPF) > tolP);
viol_Pflow(mPT) = viol_Pflow(mPT) | (dPT_pu(mPT) > tolP);

viol_Qflow(mQF) = viol_Qflow(mQF) | (dQF_pu(mQF) > tolQ);
viol_Qflow(mQT) = viol_Qflow(mQT) | (dQT_pu(mQT) > tolQ);

% Store elementwise diffs (optional, useful for debugging)
report.errcheck.diffs = struct();
report.errcheck.diffs.dPg_pu   = dPg_pu;
report.errcheck.diffs.dQg_pu   = dQg_pu;
report.errcheck.diffs.dVm_pu   = dVm;
report.errcheck.diffs.dVa_rad  = dVa_rad;
report.errcheck.diffs.dPF_pu   = dPF_pu;
report.errcheck.diffs.dPT_pu   = dPT_pu;
report.errcheck.diffs.dQF_pu   = dQF_pu;
report.errcheck.diffs.dQT_pu   = dQT_pu;

% Counts for the secondary table
report.errcheck.counts = struct( ...
    'viol_Pg', nnz(viol_Pg), ...
    'viol_Qg', nnz(viol_Qg), ...
    'viol_V',  nnz(viol_Vm), ...
    'viol_Ang',nnz(viol_Va), ...
    'viol_Pflow', nnz(viol_Pflow), ...
    'viol_Qflow', nnz(viol_Qflow));

% A table-like triplet (n_viol | max | avg) in pu/rad
report.errcheck.table = struct();
report.errcheck.table.genP  = make_all_triplet(dPg_pu,  viol_Pg);
report.errcheck.table.genQ  = make_all_triplet(dQg_pu,  viol_Qg);
report.errcheck.table.volt  = make_all_triplet(dVm,     viol_Vm);
report.errcheck.table.ang   = make_all_triplet(dVa_rad, viol_Va);

% For flows, compute severity per branch as max(PF,PT) and max(QF,QT)
Pflow_dev = nanmax([dPF_pu, dPT_pu], [], 2);
Qflow_dev = nanmax([dQF_pu, dQT_pu], [], 2);

% nanmax returns NaN if both are NaN; keep violations false there already
report.errcheck.table.Pflow = make_all_triplet(abs(Pflow_dev), viol_Pflow);
report.errcheck.table.Qflow = make_all_triplet(abs(Qflow_dev), viol_Qflow);




% ---- TRUE TRIG bus balance check (pu) ----
deltaE = Va(bi) - Va(bj) - phi_e;
Pbus = A * (Vm(bi).*Vm(bj).*absYe .* cos(deltaE));
Qbus = A * (Vm(bi).*Vm(bj).*absYe .* sin(deltaE));
Pmis = Pg_bus - Pd - Pbus;
Qmis = Qg_bus - Qd - Qbus;

viol_Pmis = abs(Pmis) > tolP;
viol_Qmis = abs(Qmis) > tolQ;
viol_bus_balance = (viol_Pmis | viol_Qmis);

% ---- TRUE TRIG branch flows and thermal/angle checks ----
r = br(on,BR_R); x = br(on,BR_X); b_total = br(on,BR_B);
rateA = br(on,RATE_A);
tapm = br(on,TAP); tapm(tapm==0)=1;
shft = br(on,SHIFT);
t = tapm .* exp(1j*deg2rad(shft));

y_series = 1./(r + 1j*x);
y_sh = 1j*(b_total);

a_f = (y_series + 0.5*y_sh)./(abs(t).^2);
b_f =  y_series ./ conj(t);
a_t = (y_series + 0.5*y_sh);
b_t =  y_series ./ t;

af_conj_r = real(conj(a_f)); af_conj_i = imag(conj(a_f));
at_conj_r = real(conj(a_t)); at_conj_i = imag(conj(a_t));
bf_mag = abs(b_f); ang_bf = angle(b_f);
bt_mag = abs(b_t); ang_bt = angle(b_t);

Vf = Vm(fbus); Vt = Vm(tbus);
d_f = Va(fbus) - Va(tbus) - ang_bf;
d_t = Va(tbus) - Va(fbus) - ang_bt;

Pf = af_conj_r.*(Vf.^2) - bf_mag.*(Vf.*Vt).*cos(d_f);
Qf = af_conj_i.*(Vf.^2) - bf_mag.*(Vf.*Vt).*sin(d_f);
Pt = at_conj_r.*(Vt.^2) - bt_mag.*(Vt.*Vf).*cos(d_t);
Qt = at_conj_i.*(Vt.^2) - bt_mag.*(Vt.*Vf).*sin(d_t);

Sf = sqrt(Pf.^2 + Qf.^2);
St = sqrt(Pt.^2 + Qt.^2);

Smax = rateA/baseMVA;
idx_rate = find(~isnan(Smax) & Smax>0);
viol_thermal_f = false(size(fbus));
viol_thermal_t = false(size(fbus));
if ~isempty(idx_rate)
    viol_thermal_f(idx_rate) = Sf(idx_rate) > (Smax(idx_rate) + tolTh);
    viol_thermal_t(idx_rate) = St(idx_rate) > (Smax(idx_rate) + tolTh);
end

% angle-difference limits (rad), tolerance tolAng
angmin = deg2rad(br(on,ANGMIN));
angmax = deg2rad(br(on,ANGMAX));
has_l = ~isnan(angmin) & (angmin > -pi);
has_u = ~isnan(angmax) & (angmax <  pi);

angdiff = Va(fbus) - Va(tbus);  % rad
viol_ang_l = false(size(fbus));
viol_ang_u = false(size(fbus));
if any(has_l), viol_ang_l(has_l) = (angdiff(has_l) < angmin(has_l) - tolAng); end
if any(has_u), viol_ang_u(has_u) = (angdiff(has_u) > angmax(has_u) + tolAng); end

% voltage & gen checks (separate tolerances)
viol_V  = (Vm < Vmin - tolV) | (Vm > Vmax + tolV);
viol_Pg = (Pg_pu_gen < Pgmin - tolPg) | (Pg_pu_gen > Pgmax + tolPg);
viol_Qg = (Qg_pu_gen < Qgmin - tolQg) | (Qg_pu_gen > Qgmax + tolQg);

% PASS/FAIL flags
pass_P  = max(abs(Pmis)) <= tolP;
pass_Q  = max(abs(Qmis)) <= tolQ;
pass_V  = ~any(viol_V);
pass_Pg = ~any(viol_Pg);
pass_Qg = ~any(viol_Qg);
pass_ang = ~(any(viol_ang_l) || any(viol_ang_u));
pass_th  = ~(any(viol_thermal_f) || any(viol_thermal_t));

report.feasibility = struct();
report.feasibility.tols = struct( ...
    'tolP_pu', tolP, 'tolQ_pu', tolQ, 'tolV_pu', tolV, ...
    'tolPg_pu', tolPg, 'tolQg_pu', tolQg, 'tolAng_rad', tolAng, 'tolTh_pu', tolTh);

report.feasibility.trig_ap = struct();
report.feasibility.trig_ap.max_abs_Pmis_pu = max(abs(Pmis));
report.feasibility.trig_ap.max_abs_Qmis_pu = max(abs(Qmis));
report.feasibility.trig_ap.pass = pass_P && pass_Q && pass_V && pass_Pg && pass_Qg && pass_ang && pass_th;
report.feasibility.trig_ap.checks = struct( ...
    'power_balance_P', pass_P, ...
    'power_balance_Q', pass_Q, ...
    'voltage_limits',  pass_V, ...
    'genP_limits',     pass_Pg, ...
    'genQ_limits',     pass_Qg, ...
    'angle_limits',    pass_ang, ...
    'thermal_limits',  pass_th);

report.feasibility.trig_ap.counts = struct( ...
    'viol_bus_balance', nnz(viol_bus_balance), ...
    'viol_voltage', nnz(viol_V), ...
    'viol_genP', nnz(viol_Pg), ...
    'viol_genQ', nnz(viol_Qg), ...
    'viol_angle', nnz(viol_ang_l | viol_ang_u), ...
    'viol_thermal', nnz(viol_thermal_f | viol_thermal_t));

% explicit violation lists (same structure as before)
report.feasibility.trig_ap.violations = struct();

bb_ids = find(viol_bus_balance);
report.feasibility.trig_ap.violations.bus_balance.bus     = bb_ids;
report.feasibility.trig_ap.violations.bus_balance.Pmis_pu = Pmis(bb_ids);
report.feasibility.trig_ap.violations.bus_balance.Qmis_pu = Qmis(bb_ids);

V_ids = find(viol_V);
report.feasibility.trig_ap.violations.voltage.bus = V_ids;

gP_ids = find(viol_Pg);
gQ_ids = find(viol_Qg);
report.feasibility.trig_ap.violations.genP.gen = gP_ids;
report.feasibility.trig_ap.violations.genQ.gen = gQ_ids;

ang_ids = find(viol_ang_l | viol_ang_u);
report.feasibility.trig_ap.violations.angle.br_on    = ang_ids;
report.feasibility.trig_ap.violations.angle.from_bus = fbus(ang_ids);
report.feasibility.trig_ap.violations.angle.to_bus   = tbus(ang_ids);

th_ids = find(viol_thermal_f | viol_thermal_t);
report.feasibility.trig_ap.violations.thermal.br_on    = th_ids;
report.feasibility.trig_ap.violations.thermal.from_bus = fbus(th_ids);
report.feasibility.trig_ap.violations.thermal.to_bus   = tbus(th_ids);

% optional top lists (same as before)
report.feasibility.trig_ap.top = build_top_lists( ...
    Pmis, Qmis, viol_V, Vm, Vmin, Vmax, ...
    viol_Pg, Pg_pu_gen, Pgmin, Pgmax, genbus, ...
    viol_Qg, Qg_pu_gen, Qgmin, Qgmax, ...
    (viol_ang_l | viol_ang_u), angdiff, angmin, angmax, fbus, tbus, ...
    (viol_thermal_f | viol_thermal_t), Sf, St, Smax, fbus, tbus, ...
    list_top);

% ---------------- Violation severities ("how bad") ----------------
viol_ang = (viol_ang_l | viol_ang_u);
viol_th  = (viol_thermal_f | viol_thermal_t);

report.feasibility.trig_ap.severity.power_balance = struct( ...
    'max_abs_Pmis_pu', max(abs(Pmis)), ...
    'avg_abs_Pmis_pu', mean(abs(Pmis)), ...
    'max_abs_Qmis_pu', max(abs(Qmis)), ...
    'avg_abs_Qmis_pu', mean(abs(Qmis)));

V_under = max(0, (Vmin - Vm));
V_over  = max(0, (Vm - Vmax));
V_margin = max(V_under, V_over);

report.feasibility.trig_ap.severity.voltage = struct( ...
    'max_margin_pu',  max(V_margin), ...
    'avg_margin_pu',  mean_if_nonempty(V_margin(viol_V)));

Pg_under = max(0, (Pgmin - Pg_pu_gen));
Pg_over  = max(0, (Pg_pu_gen - Pgmax));
Pg_margin = max(Pg_under, Pg_over);

report.feasibility.trig_ap.severity.genP = struct( ...
    'max_margin_pu',  max(Pg_margin), ...
    'avg_margin_pu',  mean_if_nonempty(Pg_margin(viol_Pg)));

Qg_under = max(0, (Qgmin - Qg_pu_gen));
Qg_over  = max(0, (Qg_pu_gen - Qgmax));
Qg_margin = max(Qg_under, Qg_over);

report.feasibility.trig_ap.severity.genQ = struct( ...
    'max_margin_pu',  max(Qg_margin), ...
    'avg_margin_pu',  mean_if_nonempty(Qg_margin(viol_Qg)));

ang_margin = zeros(size(angdiff));
if any(viol_ang_l), ang_margin(viol_ang_l) = (angmin(viol_ang_l) - angdiff(viol_ang_l)); end
if any(viol_ang_u), ang_margin(viol_ang_u) = (angdiff(viol_ang_u) - angmax(viol_ang_u)); end

report.feasibility.trig_ap.severity.angle = struct( ...
    'max_margin_rad', max(ang_margin), ...
    'avg_margin_rad', mean_if_nonempty(ang_margin(viol_ang)));

th_margin = zeros(size(Sf));
if ~isempty(idx_rate)
    exceed_f = max(0, Sf - Smax);
    exceed_t = max(0, St - Smax);
    th_margin = max(exceed_f, exceed_t);
end

report.feasibility.trig_ap.severity.thermal = struct( ...
    'max_margin_pu', max(th_margin), ...
    'avg_margin_pu', mean_if_nonempty(th_margin(viol_th)));


% ---------------- Table-style aggregate metrics (ALL elements) ----------------
report.table = struct();

% Bus balance
report.table.busP = make_all_triplet(abs(Pmis), viol_Pmis);
report.table.busQ = make_all_triplet(abs(Qmis), viol_Qmis);

% Voltage deviation from bounds (pu)
V_under = max(0, Vmin - Vm);
V_over  = max(0, Vm - Vmax);
V_dev = max(V_under, V_over);
report.table.volt = make_all_triplet(V_dev, viol_V);

% Generator P/Q deviation from bounds (pu)
Pg_under = max(0, Pgmin - Pg_pu_gen);
Pg_over  = max(0, Pg_pu_gen - Pgmax);
Pg_dev = max(Pg_under, Pg_over);
report.table.genP = make_all_triplet(Pg_dev, viol_Pg);

Qg_under = max(0, Qgmin - Qg_pu_gen);
Qg_over  = max(0, Qg_pu_gen - Qgmax);
Qg_dev = max(Qg_under, Qg_over);
report.table.genQ = make_all_triplet(Qg_dev, viol_Qg);

% Angle-difference deviation from bounds (rad)
ang_dev = zeros(size(angdiff));
if any(viol_ang_l)
    ang_dev(viol_ang_l) = abs(angmin(viol_ang_l) - angdiff(viol_ang_l));
end
if any(viol_ang_u)
    ang_dev(viol_ang_u) = abs(angdiff(viol_ang_u) - angmax(viol_ang_u));
end
viol_ang = viol_ang_l | viol_ang_u;
report.table.ang = make_all_triplet(abs(ang_dev), viol_ang);

% Thermal deviation above limits (pu)
th_dev = zeros(size(Sf));
if ~isempty(idx_rate)
    exceed_f = max(0, Sf - Smax);
    exceed_t = max(0, St - Smax);
    th_dev = max(exceed_f, exceed_t);
end
viol_th = viol_thermal_f | viol_thermal_t;
report.table.therm = make_all_triplet(abs(th_dev), viol_th);




% ---------------- PRINT (clean & unit-labeled) ----------------
print_clean_report(report, list_top);



end % main


% ================= helpers =================
function tol_out = inherit_tol(tol_in, tol_global)
if isempty(tol_in)
    tol_out = tol_global;
else
    tol_out = tol_in;
end
end

function s = stats_from_vector(v1, v2)
s = struct();
s.exp = vec_stats(v1);
s.ap  = vec_stats(v2);
if ~any(isnan([s.exp.max, s.ap.max]))
    s.diff = struct( ...
        'max_abs', max(abs(v2(:) - v1(:))), ...
        'avg_abs', mean(abs(v2(:) - v1(:))), ...
        'min_abs', min(abs(v2(:) - v1(:))) );
else
    s.diff = struct('max_abs',NaN,'avg_abs',NaN,'min_abs',NaN);
end
end

function st = vec_stats(v)
if isempty(v)
    st = struct('max',NaN,'avg',NaN,'min',NaN,'n',0);
else
    v = v(:);
    st = struct('max',max(v), 'avg',mean(v), 'min',min(v), 'n',numel(v));
end
end

function v = getfield_safe(s, path)
v = [];
try
    if numel(path)==1
        if isfield(s, path{1}), v = s.(path{1}); end
    else
        if isfield(s, path{1}) && isfield(s.(path{1}), path{2})
            v = s.(path{1}).(path{2});
        end
    end
catch
    v = [];
end
end

function cb = extract_branch_congestion(res)
cb = struct('n',0,'data',[]);
if isfield(res,'cong_branch') && isfield(res.cong_branch,'branch_idx') && ~isempty(res.cong_branch.branch_idx)
    cb.n = numel(res.cong_branch.branch_idx);
    cb.data = res.cong_branch;
end
end

function mpc_i = to_internal_copy(mpc_in)
if isfield(mpc_in, 'order')
    if ~isfield(mpc_in.order,'state') || ~strcmp(mpc_in.order.state,'i')
        mpc_i = ext2int(mpc_in);
    else
        mpc_i = mpc_in;
    end
else
    mpc_i = ext2int(mpc_in);
end
end

function top = build_top_lists(Pmis,Qmis,violV,Vm,Vmin,Vmax, ...
                               violPg,Pg,Pgmin,Pgmax,genbus, ...
                               violQg,Qg,Qgmin,Qgmax, ...
                               violAng,angdiff,angmin,angmax,fbus,tbus, ...
                               violTh,Sf,St,Smax,fbus2,tbus2, ...
                               list_top)
top = struct();

[~,ordP] = sort(abs(Pmis),'descend');
[~,ordQ] = sort(abs(Qmis),'descend');
top.bus_misP = ordP(1:min(list_top,numel(ordP)));
top.bus_misQ = ordQ(1:min(list_top,numel(ordQ)));

top.bus_V = find(violV);
if numel(top.bus_V)>list_top, top.bus_V = top.bus_V(1:list_top); end

top.genP = find(violPg); if numel(top.genP)>list_top, top.genP=top.genP(1:list_top); end
top.genQ = find(violQg); if numel(top.genQ)>list_top, top.genQ=top.genQ(1:list_top); end

top.ang = find(violAng); if numel(top.ang)>list_top, top.ang=top.ang(1:list_top); end
top.th  = find(violTh);  if numel(top.th)>list_top,  top.th=top.th(1:list_top);  end
end

function print_clean_report(report, list_top)
fprintf('\n============================================================\n');
fprintf('AC-OPF Comparison Report\n');
fprintf('============================================================\n');

% Objective
fprintf('\n[Objective] (cost units from gencost)\n');
fprintf('  f_EXP = %.6f\n', report.objective.f_exp);
fprintf('  f_AP  = %.6f\n', report.objective.f_ap);
fprintf('  df    = %.6f\n', report.objective.df);
fprintf('  gap   = %.4f %%\n', report.objective.gap_percent);

% Timing
fprintf('\n[Timing]\n');
fprintf('  Elapsed solver time [s]\n');
fprintf('    EXP: %.6f s\n', report.timing.exp.et_s);
fprintf('    AP : %.6f s\n', report.timing.ap.et_s);

if ~isnan(report.timing.exp.ipopt_iter) || ~isnan(report.timing.ap.ipopt_iter)
    fprintf('  IPOPT iterations [count] (if available)\n');
    if isnan(report.timing.exp.ipopt_iter), fprintf('    EXP: (not available)\n');
    else, fprintf('    EXP: %d\n', round(report.timing.exp.ipopt_iter)); end
    if isnan(report.timing.ap.ipopt_iter), fprintf('    AP : (not available)\n');
    else, fprintf('    AP : %d\n', round(report.timing.ap.ipopt_iter)); end
end

if ~isnan(report.timing.speedup)
    if report.timing.delta_s < 0
        fprintf('  Speed: ALL-PASS is FASTER by %.6f s (%.2f%%).  Speedup = %.3fx\n', ...
            abs(report.timing.delta_s), abs(report.timing.delta_percent), report.timing.speedup);
    elseif report.timing.delta_s > 0
        fprintf('  Speed: ALL-PASS is SLOWER by %.6f s (%.2f%%).  Speedup = %.3fx\n', ...
            report.timing.delta_s, report.timing.delta_percent, report.timing.speedup);
    else
        fprintf('  Speed: same elapsed time.  Speedup = %.3fx\n', report.timing.speedup);
    end
else
    fprintf('  Speed: not computed (missing/invalid elapsed times).\n');
end

% Errors
fprintf('\n[Variable errors] absolute difference |AP - EXP|\n');
fprintf('  Pg [MW]    : max=%9.3e | avg=%9.3e | min=%9.3e\n', report.errors.Pg.max, report.errors.Pg.avg, report.errors.Pg.min);
fprintf('  Qg [MVAr]  : max=%9.3e | avg=%9.3e | min=%9.3e\n', report.errors.Qg.max, report.errors.Qg.avg, report.errors.Qg.min);
fprintf('  Vm [pu]    : max=%9.3e | avg=%9.3e | min=%9.3e\n', report.errors.Vm.max, report.errors.Vm.avg, report.errors.Vm.min);
fprintf('  Va [deg]   : max=%9.3e | avg=%9.3e | min=%9.3e\n', report.errors.Va.max, report.errors.Va.avg, report.errors.Va.min);

fprintf('  PF [MW]    : max=%9.3e | avg=%9.3e | min=%9.3e\n', report.errors.PF.max, report.errors.PF.avg, report.errors.PF.min);
fprintf('  QF [MVAr]  : max=%9.3e | avg=%9.3e | min=%9.3e\n', report.errors.QF.max, report.errors.QF.avg, report.errors.QF.min);
fprintf('  PT [MW]    : max=%9.3e | avg=%9.3e | min=%9.3e\n', report.errors.PT.max, report.errors.PT.avg, report.errors.PT.min);
fprintf('  QT [MVAr]  : max=%9.3e | avg=%9.3e | min=%9.3e\n', report.errors.QT.max, report.errors.QT.avg, report.errors.QT.min);


% dtheta stats
fprintf('\n[Angles from solver outputs]\n');
fprintf('  dtheta_line.deg = (theta_from - theta_to) per in-service branch [deg]\n');
fprintf('    EXP: max=% .6f | avg=% .6f | min=% .6f   (n=%d)\n', ...
    report.stats.dtheta.exp.max, report.stats.dtheta.exp.avg, report.stats.dtheta.exp.min, report.stats.dtheta.exp.n);
fprintf('    AP : max=% .6f | avg=% .6f | min=% .6f   (n=%d)\n', ...
    report.stats.dtheta.ap.max, report.stats.dtheta.ap.avg, report.stats.dtheta.ap.min, report.stats.dtheta.ap.n);

fprintf('\n  deltaE_deg = (theta_i - theta_j - phi_e) over off-diagonal Ybus edges [deg]\n');
fprintf('    EXP: max=% .6f | avg=% .6f | min=% .6f   (n=%d)\n', ...
    report.stats.deltaE.exp.max, report.stats.deltaE.exp.avg, report.stats.deltaE.exp.min, report.stats.deltaE.exp.n);
fprintf('    AP : max=% .6f | avg=% .6f | min=% .6f   (n=%d)\n', ...
    report.stats.deltaE.ap.max, report.stats.deltaE.ap.avg, report.stats.deltaE.ap.min, report.stats.deltaE.ap.n);

% Congestion
fprintf('\n[Congestion from solver reports] (branches only)\n');
fprintf('  EXP: %d congested branches\n', report.congestion.exp.n);
fprintf('  AP : %d congested branches\n', report.congestion.ap.n);

if isfield(report.congestion,'match')
    cU = report.congestion.match.undirected.common_keys;
    oE = report.congestion.match.undirected.only_exp;
    oA = report.congestion.match.undirected.only_ap;

    fprintf('  Match (physical/undirected) : common=%d | only-EXP=%d | only-AP=%d\n', ...
        numel(cU), numel(oE), numel(oA));

    topK = 10;
    if ~isempty(oE)
        fprintf('  Lines congested ONLY in EXP (showing up to %d) [bus-bus]:\n', min(topK,numel(oE)));
        for i=1:min(topK,numel(oE)), fprintf('    %s\n', oE{i}); end
    end
    if ~isempty(oA)
        fprintf('  Lines congested ONLY in AP (showing up to %d) [bus-bus]:\n', min(topK,numel(oA)));
        for i=1:min(topK,numel(oA)), fprintf('    %s\n', oA{i}); end
    end
end

% Feasibility
fprintf('\n[Feasibility] AP solution evaluated under TRUE trig AC model\n');
t = report.feasibility.tols;
fprintf('  Tolerances: tolP=%.1e pu | tolQ=%.1e pu | tolV=%.1e pu | tolPg=%.1e pu | tolQg=%.1e pu | tolAng=%.1e rad | tolTh=%.1e pu\n', ...
    t.tolP_pu, t.tolQ_pu, t.tolV_pu, t.tolPg_pu, t.tolQg_pu, t.tolAng_rad, t.tolTh_pu);

fprintf('  Checks: Pbal=%s | Qbal=%s | V=%s | Pg=%s | Qg=%s | ang=%s | thermal=%s\n', ...
    yn(report.feasibility.trig_ap.checks.power_balance_P), ...
    yn(report.feasibility.trig_ap.checks.power_balance_Q), ...
    yn(report.feasibility.trig_ap.checks.voltage_limits), ...
    yn(report.feasibility.trig_ap.checks.genP_limits), ...
    yn(report.feasibility.trig_ap.checks.genQ_limits), ...
    yn(report.feasibility.trig_ap.checks.angle_limits), ...
    yn(report.feasibility.trig_ap.checks.thermal_limits));
fprintf('  PASS overall? %s\n', yn(report.feasibility.trig_ap.pass));

fprintf('  Violation counts [count]\n');
fprintf('    Bus balance (P/Q) : %d\n', report.feasibility.trig_ap.counts.viol_bus_balance);
fprintf('    Voltage           : %d\n', report.feasibility.trig_ap.counts.viol_voltage);
fprintf('    Gen P             : %d\n', report.feasibility.trig_ap.counts.viol_genP);
fprintf('    Gen Q             : %d\n', report.feasibility.trig_ap.counts.viol_genQ);
fprintf('    Angle             : %d\n', report.feasibility.trig_ap.counts.viol_angle);
fprintf('    Thermal           : %d\n', report.feasibility.trig_ap.counts.viol_thermal);

if isfield(report.feasibility.trig_ap,'violations')
    v = report.feasibility.trig_ap.violations;

    print_id_list('  Bus-balance violating buses [bus idx]:', v.bus_balance.bus, list_top);
    print_bus_balance_details(v.bus_balance, list_top);

    print_id_list('  Voltage violating buses [bus idx]:', v.voltage.bus, list_top);
    print_id_list('  GenP violating generators [gen row]:', v.genP.gen, list_top);
    print_id_list('  GenQ violating generators [gen row]:', v.genQ.gen, list_top);

    print_id_list('  Angle violating branches [on-idx]:', v.angle.br_on, list_top);
    print_branch_endpoints('  Angle violating branch endpoints [from->to]:', v.angle.from_bus, v.angle.to_bus, list_top);

    print_id_list('  Thermal violating branches [on-idx]:', v.thermal.br_on, list_top);
    print_branch_endpoints('  Thermal violating branch endpoints [from->to]:', v.thermal.from_bus, v.thermal.to_bus, list_top);
end

% ---- LaTeX-table-aligned diagnostics (ALL elements) ----
if isfield(report,'table')
    fprintf('\n[Table diagnostics: n_viol | max | avg]\n');
    fprintf('  Bus Bal. P (E-1): %4d | %.3e | %.3e\n', report.table.busP.n_viol, report.table.busP.max, report.table.busP.avg);
    fprintf('  Bus Bal. Q (E-1): %4d | %.3e | %.3e\n', report.table.busQ.n_viol, report.table.busQ.max, report.table.busQ.avg);
    fprintf('  Volt.     (E-4): %4d | %.3e | %.3e\n', report.table.volt.n_viol, report.table.volt.max, report.table.volt.avg);
    fprintf('  Gen P     (E-2): %4d | %.3e | %.3e\n', report.table.genP.n_viol, report.table.genP.max, report.table.genP.avg);
    fprintf('  Gen Q     (E-2): %4d | %.3e | %.3e\n', report.table.genQ.n_viol, report.table.genQ.max, report.table.genQ.avg);
    fprintf('  thetaij   (E-3): %4d | %.3e | %.3e\n', report.table.ang.n_viol,  report.table.ang.max,  report.table.ang.avg);
    fprintf('  Ther.     (E-2): %4d | %.3e | %.3e\n', report.table.therm.n_viol,report.table.therm.max,report.table.therm.avg);
end

% ---- Secondary table: EXP vs AP elementwise mismatch violations (pu/rad) ----
if isfield(report,'errcheck') && isfield(report.errcheck,'table')
    t2 = report.errcheck.tols;

    fprintf('\n[EXP vs AP mismatch table (pu/rad): n_viol | max | avg]\n');
    fprintf('  Tols: tolP=%.1e pu | tolQ=%.1e pu | tolV=%.1e pu | tolPg=%.1e pu | tolQg=%.1e pu | tolAng=%.1e rad\n', ...
        t2.tolP_pu, t2.tolQ_pu, t2.tolV_pu, t2.tolPg_pu, t2.tolQg_pu, t2.tolAng_rad);

    fprintf('  Gen P   (pu): %4d | %.3e | %.3e\n', report.errcheck.table.genP.n_viol,  report.errcheck.table.genP.max,  report.errcheck.table.genP.avg);
    fprintf('  Gen Q   (pu): %4d | %.3e | %.3e\n', report.errcheck.table.genQ.n_viol,  report.errcheck.table.genQ.max,  report.errcheck.table.genQ.avg);
    fprintf('  Volt    (pu): %4d | %.3e | %.3e\n', report.errcheck.table.volt.n_viol,  report.errcheck.table.volt.max,  report.errcheck.table.volt.avg);
    fprintf('  Angle  (rad): %4d | %.3e | %.3e\n', report.errcheck.table.ang.n_viol,   report.errcheck.table.ang.max,   report.errcheck.table.ang.avg);
    fprintf('  Pflow(PF/PT) : %4d | %.3e | %.3e\n', report.errcheck.table.Pflow.n_viol, report.errcheck.table.Pflow.max, report.errcheck.table.Pflow.avg);
    fprintf('  Qflow(QF/QT) : %4d | %.3e | %.3e\n', report.errcheck.table.Qflow.n_viol, report.errcheck.table.Qflow.max, report.errcheck.table.Qflow.avg);
end


fprintf('============================================================\n\n');
end

function s = yn(tf)
if tf, s='OK'; else, s='FAIL'; end
end

function v = get_scalar_or_nan(s, field)
v = NaN;
try
    if isfield(s, field) && ~isempty(s.(field))
        vv = s.(field);
        if isnumeric(vv) && isscalar(vv)
            v = double(vv);
        end
    end
catch
    v = NaN;
end
end

function print_id_list(header, ids, topk)
fprintf('%s ', header);
if isempty(ids)
    fprintf('(none)\n');
    return;
end
k = min(topk, numel(ids));
fprintf('%d shown / %d total: ', k, numel(ids));
fprintf('%d ', ids(1:k));
if numel(ids) > k, fprintf('...'); end
fprintf('\n');
end

function print_bus_balance_details(bb, topk)
if ~isfield(bb,'bus') || isempty(bb.bus)
    return;
end
k = min(topk, numel(bb.bus));
fprintf('  Bus-balance details (first %d) [Pmis/Qmis pu]:\n', k);
for i=1:k
    fprintf('    bus=%d  Pmis=%+.3e pu  Qmis=%+.3e pu\n', bb.bus(i), bb.Pmis_pu(i), bb.Qmis_pu(i));
end
end

function print_branch_endpoints(header, fb, tb, topk)
fprintf('%s ', header);
if isempty(fb) || isempty(tb)
    fprintf('(none)\n');
    return;
end
k = min(topk, numel(fb));
fprintf('%d shown / %d total:\n', k, numel(fb));
for i=1:k
    fprintf('    %d -> %d\n', fb(i), tb(i));
end
if numel(fb) > k
    fprintf('    ... and %d more\n', numel(fb)-k);
end
end

function m = mean_if_nonempty(x)
if isempty(x)
    m = NaN;
else
    m = mean(x);
end
end

function M = compare_branch_congestion(cong_exp, cong_ap)
M = struct();
M.directed  = struct('common_keys',{{}}, 'only_exp',{{}}, 'only_ap',{{}});
M.undirected= struct('common_keys',{{}}, 'only_exp',{{}}, 'only_ap',{{}});

if isempty(cong_exp) || ~isfield(cong_exp,'n') || cong_exp.n==0 || isempty(cong_exp.data) || ...
   isempty(cong_ap)  || ~isfield(cong_ap,'n')  || cong_ap.n==0  || isempty(cong_ap.data)
    return;
end

cbE = cong_exp.data;
cbA = cong_ap.data;

req = {'from_bus','to_bus'};
for k=1:numel(req)
    if ~isfield(cbE,req{k}) || ~isfield(cbA,req{k})
        return;
    end
end

keysE_dir = arrayfun(@(i) sprintf('%d-%d', cbE.from_bus(i), cbE.to_bus(i)), (1:numel(cbE.from_bus))', 'uni', false);
keysA_dir = arrayfun(@(i) sprintf('%d-%d', cbA.from_bus(i), cbA.to_bus(i)), (1:numel(cbA.from_bus))', 'uni', false);

keysE_und = arrayfun(@(i) sprintf('%d-%d', min(cbE.from_bus(i),cbE.to_bus(i)), max(cbE.from_bus(i),cbE.to_bus(i))), ...
    (1:numel(cbE.from_bus))', 'uni', false);
keysA_und = arrayfun(@(i) sprintf('%d-%d', min(cbA.from_bus(i),cbA.to_bus(i)), max(cbA.from_bus(i),cbA.to_bus(i))), ...
    (1:numel(cbA.from_bus))', 'uni', false);

common_dir = intersect(keysE_dir, keysA_dir, 'stable');
onlyE_dir  = setdiff(keysE_dir, keysA_dir, 'stable');
onlyA_dir  = setdiff(keysA_dir, keysE_dir, 'stable');

common_und = intersect(keysE_und, keysA_und, 'stable');
onlyE_und  = setdiff(keysE_und, keysA_und, 'stable');
onlyA_und  = setdiff(keysA_und, keysE_und, 'stable');

M.directed.common_keys   = common_dir;
M.directed.only_exp      = onlyE_dir;
M.directed.only_ap       = onlyA_dir;

M.undirected.common_keys = common_und;
M.undirected.only_exp    = onlyE_und;
M.undirected.only_ap     = onlyA_und;
end

function T = make_all_triplet(x, viol)
% count of violations, max and avg over ALL elements
x = x(:);
viol = viol(:);

T = struct();
T.n_viol = nnz(viol);

if isempty(x)
    T.max = NaN;
    T.avg = NaN;
else
    T.max = max(x);
    T.avg = mean(x);
end
end






