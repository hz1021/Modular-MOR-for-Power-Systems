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

- Simulink, and Simulink Control Design for model linearisation.
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
        ↓
plotTest
```

`main_sub_118_bus` prepares the nonlinear 118-bus operating point and model data required by the Simulink validation.

`plotFun` runs the time-domain validation examples and generates the corresponding comparison figures.

`plotTest` runs the frequency-domain validation examples and generates the computational comparison results.

## Project files

- `modularMOR_main.m` — main three-subsystem modular MOR and evaluation script.
- `IBR_modelGenerator.m` — collects the external-area subsystems, and constructs the dynamic interface and interconnected models.
- `DK_iteration_adapt.m` — robust-performance calculation used to derive local reduction accuracy requirements for subsystems.
- `greedyLoewner.m` — greedy Loewner method.
- `mimoLoewnerWeightedAAA_stableStrict.m` — Loewner framework implementation, which produces stable ROM and applies to MIMO cases.
- `linSys1.mat`, `linSys2.mat`, and `linSys3.mat` — full-order linearised subsystem models at the investigated operating mode.
- `main_sub_118_bus.m` — initialises the nonlinear 118-bus case and linearises the Simulink model.
- `subsystem_118_bus.slx` — main nonlinear 118-bus Simulink model.
- `subsystem_118_bus_self_clearing_blocking.slx` — Simulink model for the self-clearing inverter-blocking example.
- `plotFun.m` — runs the time-domain validation cases and generates the figures.
- `case118_system*.xlsx` — 118-bus case-study network data.
- `nominal*.mat` — equilibrium data used by the Simulink validation.
- `startup_project.m` — MATLAB path and Simulink cache/code-generation setup.

The generated files

```text
sys_j_BT.mat
sys_j_AA.mat
```

contain the subsystem ROMs obtained from the (frequency weighted) balanced truncation and greedy Loewner workflows, respectively.

## 118-bus validation

The nonlinear 118-bus case study is initialized using

```matlab
main_sub_118_bus
```

and the time-domain/frequency-domain simulations are run using

```matlab
plotFun
plotTest
```

The provided validation cases include:

- an active-power reference perturbation;
- a self-clearing grid-following inverter blocking event.

For the self-clearing case, grid-following inverter 2 in subsystem 3 is blocked at 0.10 s and restored at 0.11 s, corresponding to a 10 ms blocking interval. The event is implemented as temporary suppression of the selected converter injection followed by recovery, rather than as a permanent topology change.

## Reusing generated ROMs

After `modularMOR_main.m` has been run once, `sys_j_BT.mat` and `sys_j_AA.mat` can be retained and reused when only the Simulink validation needs to be repeated. This avoids repeating the full robust-performance and MOR computation.

## Maintainer & Contributing

- Maintainer: **Hanqing Zhang** - [MAC-X Lab](https://giordanoscarciotti.com/mac-x-lab/) and **Pudong Ge** - [Control and Power Group](https://profiles.imperial.ac.uk/pudong.ge19).
- Contributions and bug reports are welcome via issues or pull requests.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Citation

If you use this repository in academic work, please cite the associated publication. The full citation will be added once the final bibliographic information is available.

## Acknowledgments

The authors would like to thank [Dr. Mohammad Fahim Shakib](https://www.tue.nl/en/research/researchers/fahim-shakib), [Dr. Luuk Poort](https://www.linkedin.com/in/luuk-poort/), and [Dr. Lars A. L. Janssen](https://www.linkedin.com/in/lars-janssen/) for their dedicated help with building up the modular MOR framework given in this project.
