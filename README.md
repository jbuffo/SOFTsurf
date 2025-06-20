# SOFTsurf
This is a code to automatically generate lookup functions for ice salinity (S_ice) as a function of ocean salinity (S_oc) and ice-ocean interface thermal gradient (dT/dz) for binary salt solutions (e.g., NaCl, MgSO4, etc.).

The intention of this model is for use as an efficient parameterization of salt entrainment at ice-ocean interfaces, particularly for applications to ice-ocean worlds (e.g., Europa, Enceladus, etc.).

This readme will describe the prerequisites needed to run the model, the basic functionality of the model (how to run it), the model's outputs, and how to implement them.

## Prerequisites
This model uses the existing model SOFTBALL (SOLidification, Flow and Thermodynamics in Binary ALLoys) described in [Parkinson et al., (2020)](https://www.sciencedirect.com/science/article/pii/S2590055219300599). Instructions for downloading SOFTBALL and its supporting software (Chombo) can be found [HERE](https://github.com/jrgparkinson/mushy-layer/tree/master).

The model also uses MATLAB, including its Parallel Computing Toolbox.

## Downloading SOFTsurf
```bash
$ git clone https://github.com/jbuffo/SOFTsurf.git SOFTsurf
```
## Running SOFTsurf
SOFTsurf can be run by editing the run_SOFTsurf.m MATLAB script and then either running that script from a terminal (reccomended) or running in MATLAB (less ideal for cluster usage)

OR

Can be called directly from a terminal as a function.

In either case, SOFTsurf is a heavily parallelized code that uses the 'parfor' capabilities of Matlab. The model will distribute instances of SOFTBALL simulations across multiple cores (up to the number of 'salinity' values defined).

### Editing run_SOFTsurf.m and running in the terminal
The Matlab script simply defines the following variables needed by the code and then runs the script. Edit the variables for your specific binary salt solution.

| Variable   | Description                             | MATLAB Type        | Example Value                                      |
|------------|-----------------------------------------|--------------------|----------------------------------------------------|
| `salinity` | Salinity values to simulate (ppt)       | array (1×N double) | `[10, 25, 50, 75, 100, 125, 150, 175, 200, 225]`   |
| `slope`    | Liquidus slope for salt solution (K/ppt)| scalar (double)    | `-0.0913`                                          |
| `C_e`      | Eutectic concentration (ppt)            | scalar (double)    | `230`                                              |
| `beta`     | Haline contraction coefficient (1/ppt)  | scalar (double)    | `7.6e-4`                                           |
| `g`        | Gravity (m/s²)                          | scalar (double)    | `1.32`                                             |
| `k_s`      | Salt diffusivity in water (m²/s)        | scalar (double)    | `1.5e-9`                                           |
| `saltname` | Name of salt (used in filenames)        | character array    | `'NaCl'`                                           |
| `path`     | Path to SOFTBALL download               | character array    | `'/path/to/file'`                                  |

Then in terminal:
```bash
$ matlab -batch "run_SOFTsurf"
```

### Running directly from the terminal
```bash
$ matlab -batch "makeSurface(salinity,slope,C_e,beta,g,k_s,saltname,path);"
```

## Model outputs and application
Primary model outputs include:
- A file named 'All_values_array.mat' containing the raw data in a m x n loadable cell array where each row (m) has n different resolution simulations for each salinity value. m is the number of salinities specified. Currently n=4, and are predefined resolution simulations, more could be added as desired to the makeSurface.m function if more resolution/coverage is desired.
- A saved figure to check the goodness of fit of the surface to the simulated data.
- A file named 'Surf(saltname).mat' containing the input parameters and a function object called 'SOFTsurf'.

The 'SOFTsurf' object includes a function handle 'S_ice', which can be used to acquire interpolated ice salinity values from the generated S_oc vs. dT/dz vs. S_ice surface.

This is done using the syntax: SOFTsurf.S_ice(S_oc, dT_dz), where:
  
   S_oc   = ocean salinity value
   
   dT_dz  = thermal gradient at the freezing front


