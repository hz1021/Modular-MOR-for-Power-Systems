# Modular Model Order Reduction for Power Systems

MATLAB scripts for reducing interconnected power-system models while preserving their frequency-domain behavior. The primary workflow splits a three-subsystem model into two external areas and a study area, computes local accuracy bounds from a global robust-performance specification, reduces the subsystems, and reconnects them. It includes balanced truncation and a greedy MIMO Loewner/AAA reduction implementation.

## Why use this project

- Applies modular reduction to interconnected systems instead of reducing only the assembled model.
- Translates a global frequency-weighted error requirement into subsystem-level specifications.
- Includes balanced truncation and a stable, frequency-weighted Loewner/AAA method for comparison.
- Includes a separate 118-bus Simulink and power-flow case-study workflow.

## Requirements

The main reduction workflow requires:

- MATLAB with Control System Toolbox and Robust Control Toolbox.
- Parallel Computing Toolbox, because the D-K bound calculation uses `parfor` and starts a parallel pool.
- [YALMIP](https://yalmip.github.io/) and a working [MOSEK](https://www.mosek.com/) installation/license. `DK_iteration_adapt.m` configures YALMIP to use MOSEK.

The optional 118-bus workflow additionally requires Simulink, Simulink Control Design (for `linearize`), and [MATPOWER](https://matpower.org/) on the MATLAB path. MATLAB and third-party toolbox versions are not specified in this repository; use releases compatible with the functions above.

## Get started

1. Clone or download the repository and open MATLAB in its root directory.
2. Install and configure the required MATLAB products, YALMIP, and MOSEK. Confirm that `yalmip` and `mosek` are available to MATLAB.
3. Run the main reduction workflow:

   ```matlab
   startup_project
   modularMOR_main
   ```

   `startup_project` adds the project folders to the MATLAB path and places Simulink cache/code-generation files in `simulink_files/`. The reduction script uses the bundled `linSys1.mat`, `linSys2.mat`, and `linSys3.mat` models. It runs the robust-performance bound calculation, balanced-truncation and Loewner/AAA reductions, and plots the interconnected model and weighted error. It also saves `sys_j_BT.mat` and `sys_j_AA.mat` in the current working directory.

The optional 118-bus scripts are `main_sub_118_bus.m` and `plotFun.m`. They load the supplied case workbooks and Simulink models; run `startup_project` first, and ensure MATPOWER and the required Simulink products are installed. These workflows may run a power flow, linearize or simulate a model, and use additional model-specific data files.

## Project files

- `modularMOR_main.m`: main three-subsystem model reduction and evaluation script.
- `IBR_modelGenerator.m`, `DK_iteration_adapt.m`, and `greedyLoewner.m`: interconnection construction, robust-performance bound calculation, and Loewner reduction workflow.
- `mimoLoewnerWeightedAAA_stableStrict.m`: stable, frequency-weighted MIMO Loewner/AAA implementation.
- `linSys*.mat`: bundled linear subsystem models used by the main workflow.
- `main_sub_118_bus.m`, `plotFun.m`, and `subsystem_118_bus*.slx`: 118-bus case-study scripts and Simulink models.
- `case118_system*.xlsx` and `nominal*.mat`: case-study input and equilibrium data.
- `startup_project.m`: MATLAB path and Simulink cache setup.

## Help and contributing

For project questions or bug reports, open an issue in the repository where you obtained this code and include the MATLAB release, installed toolbox versions, solver details, the command you ran, and the complete error message. Contributions can be proposed as pull requests; please describe the affected workflow and include reproducible validation steps. No separate contribution guide is included.

The license is [MIT](LICENSE). The license file credits Hanqing Zhang; no separate maintainer contact information is provided in this repository.