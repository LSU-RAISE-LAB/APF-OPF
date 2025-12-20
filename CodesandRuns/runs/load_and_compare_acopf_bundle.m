function report = load_and_compare_acopf_bundle(matfile, compare_mode, tol)
%LOAD_AND_COMPARE_ACOPF_BUNDLE Load saved ACOPF bundle and run comparison
%
% report = load_and_compare_acopf_bundle(matfile, compare_mode, tol)
%
% Inputs
%   matfile       : path to .mat file produced by run_acopf_bundle_save
%   compare_mode  : string specifying comparison:
%                   'exp_vs_ap_dcpf'
%                   'exp_vs_ap_dcopf'
%                   'ap_dcpf_vs_ap_dcopf'
%   tol           : (optional)
%                   - scalar (legacy behavior, same tol for all)
%                   - struct with fields:
%                       tol.PQ, tol.V, tol.Pg, tol.Qg,
%                       tol.angle, tol.thermal
%
% Output
%   report        : struct returned by compare_acopf_results

    % ---------------- Load bundle ----------------
    assert(exist(matfile,'file')==2, 'MAT file not found.');
    S = load(matfile);
    assert(isfield(S,'out'), 'MAT file does not contain struct "out".');

    out = S.out;

    % ---------------- Select comparison ----------------
    switch lower(compare_mode)
        case 'exp_vs_ap_dcpf'
            res1 = out.exp;       ok1 = out.ok.exp;
            res2 = out.ap_dcpf;   ok2 = out.ok.ap_dcpf;
            label = 'EXP vs ALL-PASS (DC-PF)';

        case 'exp_vs_ap_dcopf'
            res1 = out.exp;       ok1 = out.ok.exp;
            res2 = out.ap_dcopf;  ok2 = out.ok.ap_dcopf;
            label = 'EXP vs ALL-PASS (DC-OPF)';

        case 'ap_dcpf_vs_ap_dcopf'
            res1 = out.ap_dcpf;   ok1 = out.ok.ap_dcpf;
            res2 = out.ap_dcopf;  ok2 = out.ok.ap_dcopf;
            label = 'ALL-PASS (DC-PF) vs ALL-PASS (DC-OPF)';

        otherwise
            error('Invalid compare_mode.');
    end

    fprintf('\n============================================================\n');
    fprintf('Loaded file : %s\n', matfile);
    fprintf('Case        : %s\n', out.meta.case_name);
    fprintf('Comparison  : %s\n', label);
    fprintf('============================================================\n');

    % ---------------- Handle tolerances ----------------
    if nargin < 3 || isempty(tol)
        % default (safe & tight)
        tol = struct();
        tol.PQ      = 1e-3;
        tol.V       = 1e-4;
        tol.Pg      = 1e-4;
        tol.Qg      = 1e-4;
        tol.angle   = 1e-4;
        tol.thermal = 1e-3;
    end



    % Map tol struct -> compare_acopf_results name-value pairs
report = compare_acopf_results(res1, ok1, res2, ok2, ...
    'tolP',   tol.PQ, ...
    'tolQ',   tol.PQ, ...
    'tolV',   tol.V, ...
    'tolPg',  tol.Pg, ...
    'tolQg',  tol.Qg, ...
    'tolAng', tol.angle, ...
    'tolTh',  tol.thermal);


    % Attach tol struct explicitly for traceability
    report.meta = struct();
    report.meta.case_name = out.meta.case_name;
    report.meta.compare_mode = compare_mode;
    report.meta.tol = tol;

end
