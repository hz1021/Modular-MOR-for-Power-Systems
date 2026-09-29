# Modular Model Order Reduction for Power Systems

MATLAB/Simulink scripts for reducing interconnected power system models while preserving their dynamical behaviour. The primary workflow considers a three-subsystem model consisting of two external areas and one study area, computes local reduction accuracy requirements from a global reduction accuracy specification, reduces the subsystems independently, and reconnects them.

The repository includes two reduction approaches: (frequency weighted) balanced truncation and a greedy Loewner method. It also contains a nonlinear 118-bus Simulink case study for time-domain validation.

## Why use this project

- Applies modular reduction to interconnected systems rather than reducing only the assembled model.
- Translates a global frequency-weighted approximation requirement into subsystem-level specifications.
- Includes (frequency-weighted) balanced truncation and a greedy Loewner method for comparison.
- Reassembles the subsystem ROMs and evaluates the resulting interconnected approximation.
- Includes a separate nonlinear 118-bus Simulink and power-flow validation workflow.
- Includes both reference-perturbation and self-clearing inverter-blocking time-domain examples.

## Requirements

The main modular MOR workflow requires:

- MATLAB with Control System Toolbox and Robust Control Toolbox.
- Parallel Computing Toolbox, because the robust-performance calculation uses `parfor`.
- [YALMIP](https://yalmip.github.io/) and a working [MOSEK](https://www.mosek.com/) installation/license.

The 118-bus validation workflow additionally requires:

- Simulink;
- Simulink Control Design, for model linearisation;
- [MATPOWER](https://matpower.org/) on the MATLAB path.

The supplied Simulink models were saved using MATLAB/Simulink R2025b.

## Get started

1. Clone or download the repository and open MATLAB in the its root directory.
   
2. Install and configure the required MATLAB products, YALMIP, and MOSEK. Confirm that `yalmip` and `mosek` are available to MATLAB.

3. Run the main modular model order reduction workflow

```matlab
modularMOR_main
```
`modularMOR_main` performs the robust-performance calculation, computes the local subsystem requirements, applies the two MOR approaches, reconnects the reduced subsystem models, and evaluates the resulting interconnected ROMs.

It also generates

```text
sys_j_BT.mat
sys_j_AA.mat
```

in the current working directory. These are generated outputs of the MOR workflow, not input files required for the first run.

4. For a complete first-time reproduction, run the following scripts sequentially in the **same MATLAB session**:

```matlab
main_sub_118_bus
plotFun
plotTest
```

The workflow is therefore

```text
modularMOR_main
        ↓
generate sys_j_BT.mat and sys_j_AA.mat
        ↓
main_sub_118_bus
        ↓
plotFun
&darr;
plotTest
```

`main_sub_118_bus` prepares the nonlinear 118-bus operating point and model data required by the Simulink validation.

`plotFun` runs the time-domain validation examples and generates the corresponding comparison figures.

`plotTest` runs the frequency-domain validation examples and generates the computational comparison results.

## Project files

- `modularMOR_main.m` — main three-subsystem modular MOR and evaluation script.
- `IBR_modelGenerator.m` — constructs the subsystem, interconnection, and interconnected models.
- `DK_iteration_adapt.m` — robust-performance calculation used to derive local subsystem accuracy requirements.
- `greedyLoewner.m` — greedy point-selection procedure used in the Loewner reduction workflow.
- `mimoLoewnerWeightedAAA_stableStrict.m` — stable, frequency-weighted MIMO Loewner/AAA implementation.
- `linSys1.mat`, `linSys2.mat`, and `linSys3.mat` — full-order linear subsystem models used as inputs to the modular MOR workflow.
- `main_sub_118_bus.m` — initializes the nonlinear 118-bus case and linearises the Simulink model.
- `subsystem_118_bus.slx` — main nonlinear 118-bus Simulink model.
- `subsystem_118_bus_self_clearing_blocking.slx` — Simulink model for the temporary inverter-blocking example.
- `plotFun.m` — runs the time-domain validation cases and generates the figures.
- `case118_system*.xlsx` — 118-bus case-study network data.
- `nominal*.mat` — equilibrium data used by the Simulink validation.
- `startup_project.m` — MATLAB path and Simulink cache/code-generation setup.

The generated files

```text
sys_j_BT.mat
sys_j_AA.mat
```

contain the subsystem ROMs obtained from the balanced-truncation and Loewner/AAA workflows, respectively.

## 118-bus validation

The nonlinear 118-bus case study is initialized using

```matlab
main_sub_118_bus
```

and the time-domain simulations are run using

```matlab
plotFun
```

The provided validation cases include:

- an active-power reference perturbation;
- a self-clearing grid-following inverter blocking event.

For the self-clearing case, grid-following inverter 2 in subsystem 3 is blocked at `0.10 s` and restored at `0.11 s`, corresponding to a 10 ms blocking interval. The event is implemented as temporary suppression of the selected converter injection followed by recovery, rather than as a permanent topology change.

## Reusing generated ROMs

After `modularMOR_main.m` has been run once, `sys_j_BT.mat` and `sys_j_AA.mat` can be retained and reused when only the Simulink validation needs to be repeated.

For example:

```matlab
startup_project
config_file

load linSys1.mat
load linSys2.mat
load linSys3.mat

G = {linSys1; linSys2; linSys3};
[sys_j, ~, ~] = IBR_modelGenerator(G);

load sys_j_BT.mat
load sys_j_AA.mat

main_sub_118_bus
plotFun
```

This avoids repeating the full robust-performance and MOR computation.

## Maintainer & Contributing

- Maintainer: **Hanqing Zhang** — [MAC-X Lab](https://giordanoscarciotti.com/mac-x-lab/).
- Contributions and bug reports are welcome via issues or pull requests.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Citation

If you use this repository in academic work, please cite the associated publication. The full citation will be added once the final bibliographic information is available.
