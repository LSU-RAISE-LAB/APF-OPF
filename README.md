APF-OPF: A Solver-Oriented Fractional Framework for AC Optimal Power Flow






Official code, datasets, and numerical results for a solver-oriented AC Optimal Power Flow framework.

🔬 Overview

This repository provides implementations and numerical evaluation tools for an alternative AC Optimal Power Flow (AC OPF) formulation developed to improve numerical behavior when solved with interior-point nonlinear programming methods.

The framework is designed to integrate seamlessly with standard OPF toolchains and enables systematic evaluation of solver performance, robustness, and scalability across a wide range of benchmark power-system networks.

🌟 Features

Solver-Compatible AC OPF Models
Formulations designed for use with interior-point solvers

Multiple Modeling Environments
Support for both YALMIP and CasADi backends

Scalable Evaluation
Tested on IEEE, MATPOWER, and PGLib transmission networks, including large-scale systems

Unified Experiment Pipeline
Automated routines for:

runtime and iteration comparison

feasibility verification

voltage, angle, and congestion consistency analysis

Reproducible Results
All numerical results can be regenerated using the provided scripts and case lists

📚 Citation

If you use the content of this repository in your research, please cite the corresponding work when available.

Citation information will be added upon publication.

👤 Author

Milad Hasanzadeh
Department of Electrical and Computer Engineering
Louisiana State University

📧 mhasa42@lsu.edu

📅 Release Date: 2025
📄 License: Academic and Research Use Only
