function report = run_acopf_compare(case_name, mode, zar, loading_condition_per)

%RUN_ACOPF_COMPARE Run EXP vs ALL-PASS AC-OPF and compare results
% Automatically saves Command Window output to a .txt file

    % ----------- Input defaults & sanity checks -----------
    if nargin < 2
        error('Usage: run_acopf_compare(case_name, mode, [zar])');
    end
    if nargin < 3 || isempty(zar)
        zar = 0;     % default: no RATE_A scaling
    end

    if nargin < 4 || isempty(loading_condition_per)
        loading_condition_per = 100;   % default: no load scaling
    end


    % ----------- Prepare log directory & file -----------
    if ischar(case_name) || isstring(case_name)
        cname = char(case_name);
    else
        cname = 'mpc';
    end

    logdir = 'logs';
    if ~exist(logdir,'dir')
        mkdir(logdir);
    end

    timestamp = datestr(now,'yyyy-mm-dd_HHMMSS');
    logfile = fullfile(logdir, ...
        sprintf('ACOPF_%s_mode%d_zar%.2f_%s.txt', cname, mode, zar, timestamp));

    diary(logfile);
    diary on;

    fprintf('\n=============================================\n');
    fprintf(' AC-OPF COMPARISON LOG START\n');
    fprintf(' Case       : %s\n', cname);
    fprintf(' Mode       : %d\n', mode);
    fprintf(' RATE_A zar : %.2f %%\n', zar);
    fprintf(' Timestamp  : %s\n', timestamp);
    fprintf(' Log file   : %s\n', logfile);
    fprintf(' Load level : %.1f %%\n', loading_condition_per);

    fprintf('=============================================\n\n');

    try
        % ----------- Parameters -----------
        a_allpass  = 0.5;
        ipopt_opts = struct();
        dctype     = 0;

        % ----------- Run OPF variants -----------
        switch mode
            case 1
                % ===== CasADi + IPOPT =====
                fprintf('> [CasADi+IPOPT] Running EXP (admittance-angle included)...\n');
                [res1, ok1] = acopf_exp_ipopt_casadi( ...
                    case_name, zar, ipopt_opts, loading_condition_per);


                fprintf('> [CasADi+IPOPT] Running ALL-PASS (admittance-angle included)...\n');
                [res2, ok2] = acopf_allpass_ipopt_casadi( ...
                    case_name, zar, a_allpass, ipopt_opts, dctype, [], loading_condition_per);


            case 2
                % ===== YALMIP + IPOPT =====
                fprintf('> [YALMIP+IPOPT] Running EXP (admittance-angle included)...\n');
                [res1, ok1] = acopf_exp( ...
                    case_name, zar, ipopt_opts, loading_condition_per);


                fprintf('> [YALMIP+IPOPT] Running ALL-PASS (admittance-angle included)...\n');
                [res2, ok2] = acopf_allpass( ...
                    case_name, zar, a_allpass, ipopt_opts, dctype, loading_condition_per);


            otherwise
                error('Invalid mode. Use 1 (CasADi+IPOPT) or 2 (YALMIP+IPOPT).');
        end

        % ----------- Compare results -----------
        fprintf('\n> Comparing results (EXP vs ALL-PASS)...\n');

        report = compare_acopf_results(res1, ok1, res2, ok2, ...
            'tolP', 1e-1, ...
            'tolQ', 1e-1, ...
            'tolV', 1e-4, ...
            'tolPg', 1e-4, ...
            'tolQg', 1e-4, ...
            'tolAng', 1e-4, ...
            'tolTh', 1e-2);


        fprintf('\n=============================================\n');
        fprintf(' COMPLETED SUCCESSFULLY\n');
        fprintf('=============================================\n\n');

    catch ME
        fprintf('\n=============================================\n');
        fprintf(' ERROR OCCURRED\n');
        fprintf(' Message: %s\n', ME.message);
        fprintf('=============================================\n\n');

        diary off;
        rethrow(ME);
    end

    diary off;
end
