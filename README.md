# Steady Flow Simulation

This project simulates steady, incompressible flow over a cambered airfoil using a simple panel method based on singularity distributions. It calculates the velocity distribution, pressure coefficient, total lift, and moment for a given airfoil geometry and flow conditions.

## Features

- Custom airfoil geometry generation with camber and thickness
- Discretization into singularity panels
- Influence coefficient matrix calculation
- Solving linear system for circulation strengths
- Calculation of lift and moment coefficients
- Tabulated output for clarity

## Airfoil Configuration

- Camber: 2% of chord
- Maximum camber position: 40% of chord
- Thickness: 12% of chord
- Chord length: 1.0 m
- Panels: 8 (discretization points)
