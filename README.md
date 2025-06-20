# SOFTsurf
This is a code to automatically generate lookup functions for ice salinity (S_ice) as a function of ocean salinity (S_oc) and ice-ocean interface thermal gradient (dT/dz) for binary salt solutions (e.g., NaCl, MgSO4, etc.).

The intention of this model is for use as an efficient parameterization of salt entrainment at ice-ocean interfaces, particularly for applications to ice-ocean worlds (e.g., Europa, Enceladus, etc.).

This readme will describe the prerequisites needed to run the model, the basic functionality of the model (how to run it), the model's outputs, and how to implement them.

## Prerequisites
This model uses the existing model SOFTBALL (SOLidification, Flow and Thermodynamics in Binary ALLoys) described in [Parkinson et al., (2020)](https://www.sciencedirect.com/science/article/pii/S2590055219300599). Instructions for downloading SOFTBALL and its supporting software (Chombo) can be found [HERE](https://github.com/jrgparkinson/mushy-layer/tree/master).

## Downloading SOFTsurf
```bash
$ git clone https://github.com/jbuffo/SOFTsurf.git SOFTsurf
