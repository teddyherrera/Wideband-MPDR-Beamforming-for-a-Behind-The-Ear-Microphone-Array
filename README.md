# Wideband MPDR Beamforming Simulation for Hearing Aid Microphone Arrays
This repository contains a MATLAB simulation framework for FIR-based  Minimum Power Distortionless Response (MPDR) beamforming applied to a behind-the-ear (BTE) hearing aid microphone array.

The project models space–time beamforming using finite impulse response (FIR) filters and evaluates array performance under different constraint configurations (look direction only, optional null constraints).

The simulation follows a space–time formulation of wideband beamforming, consistent with classical treatments in array signal processing, and is intended for performance analysis, visualization, and experimentation, rather than real-time deployment.

## Simulation Overview
1. Define physical and signal parameters
    * Speed of Sound
    * Sampling frequency
    * FIR filter length
    * Signal duration
2. Define Microphone Array Geometry
    * Microphone positions in 3D
    * Inter-microphone spacing representative of BTE devices
3. Steering Vector Construction
    * Space-time steering vector for the desired look direction
    * Optional space-time steering vectors for null directions
4. Constraint Matrix Assembly
    * Distortionless response constraint for target direction
    * Optional null constraints for interference suppression
5. MPDR Weight Computation
    * FIR beamformer weights computed using space-time covariance modeling
    * Supports zero or multiple nulls
6. Beam Pattern Evaluation
    * Directional response visualization
    * Used to assess spatial selectiity and null placement capability

## How to Run
1. Open MATLAB 2021b or later
2. Run:

        hearingAidSim
3. Inspect generated beam patterns and weight structures

References
* H. L. Van Trees, Optimum Array Processing, Wiley
* Kayser et al., “Head-related transfer functions of hearing aid microphones,” 2009
