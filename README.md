# SOFTsurf
This is a code to automatically generate lookup functions for ice salinity (S_ice) as a function of ocean salinity (S_oc) and ice-ocean interface thermal gradient (dT/dz) for binary salt solutions (e.g., NaCl, MgSO4, etc.).

The intention of this model is for use as an efficient parameterization of salt entrainment at ice-ocean interfaces, particularly for applications to ice-ocean worlds (e.g., Europa, Enceladus, etc.).

This readme will describe the prerequisites needed to run the model, the basic functionality of the model (how to run it), the model's outputs, and how to implement them.

## Prerequisites
This model uses the existing model SOFTBALL (SOLidification, Flow and Thermodynamics in Binary ALLoys) described in [Parkinson et al., (2020)](https://www.sciencedirect.com/science/article/pii/S2590055219300599). Instructions for downloading SOFTBALL and its supporting software (Chombo) can be found [HERE](https://github.com/jrgparkinson/mushy-layer/tree/master).

The model also uses Matlab.

## Downloading SOFTsurf
```bash
$ git clone https://github.com/jbuffo/SOFTsurf.git SOFTsurf
```
## Running SOFTsurf
SOFTsurf can be run by editing the run_SOFTsurf.m Matlab script and then either running that script from a terminal (reccomended) or running in Matlab

OR

can be called directly from a terminal as a function.

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


