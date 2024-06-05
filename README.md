# AvianDynStab

This repository includes all of the necessary codes used in the article: TBD.

## Citation

## Abstract

## Introduction

## Overview
run script.py -> calculate the longituidnal stability of a rigid configuration
    - relies on a source of aerodynamic coefficient data

aero_functions.py
    - finds trim state based on conditions and morphing parameters and aerodynamic data
        - exits early if no trim state solution is found
        - optimizer solves for angle of attack, airspeed, and flight path angle
        - extracts aerodynamic coefficients from input configuration and data
analyze_data.r
    - relies on "function_data.RData"
    - relies on "LongDynStability_Rigid.csv", which comes from run_script.py. This file contains the flight condition and parameter values with the eigenvalues of the system.
        - removes nonphysical configurations (line 64)
        - filters the data to only consider certain dihedral and forward sweep angles (line 64)
        - extracts additional configuration information from the filtered data (line 83)
        - find the short period and phugoid mode relevant configurations (line 119)
        - creates linear models as functions of morphing inputs (line 124)
        - relabels data per changing manus angle (??)
    - does a lot, compares prediction data against numerical and experimental data

subsample_gullwings.R
    - relies on oriented wing information
        - limit one sample per degree
        - creates a subsample of the gullwing data that houses a row per single degree of deflection
        - saves that data to a new file
    - relies on the subsampled file from above
        - subsampled frames for the wing (and tail?) 

determine_functions.R
    - relies on aerodynamic data
    - runs data based on shoulder joints for a specific pitch rate
    Inertia Data
       - adjusts the model inertia per shoulder position
       - models changes in x and z cg as a function of shoulder position for each elbow and manus angle
       - 

extract_stats.R
    - isolates the effects of each morphing parameter on the stability modes
        - effect of wrist
        - effect of elbow
        - effect of sweep
    - defines static stability about the trim state (find number of stable and unstable configurations)
    - defines dynamics stabiity about the trim state (find number of stable and unstable configurations)
    - extracts flying qualities approximations
    - extracts gust response characteristics (both stable and unstable)


## Flow


