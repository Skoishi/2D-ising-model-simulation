# 2D-ising-model-simulation

# Disorder Effect in 2D Ising Model

## Overview

This project explores the disorder effect in a 2D Ising model, implemented in C++. The main code and user interface are provided in `2d_ising.cpp`, compiled using `Make_2d_ising`. The project includes classes for calculating and storing configuration status with and without disorder.

### Main Components

- **Magnet class & Matrix class:**
  - `Magnet.h`
  - `Magnet.cpp`

  These classes handle configurations without disorder, providing the foundation for the calculations.

- **Disorder class:**
  - `Disorder.h`
  - `Disorder.cpp`

  This class deals with configurations including disorder and explores the impact of disorder on the Ising model.

- **Plotting class:**
  - `2d_ising_visual.ipynb` (for plotting configurations)
  - `Energy_calc.ipynb` (for plotting energy)

### Results

- `E_overflip.jpg`: Visualization showing the impact of disorder on the system.
- `DE_vs_C.jpg`: Graph displaying the relationship between energy and concentration.
- `DE_vs_C_log.jpg`: Logarithmic version of the previous graph.

## Physics and Data Analysis

### Animation Observations

Observations from the animation include:

- Introduction of disorder enables additional local minimum energy states.
- At low temperatures, even with 40% disorder, the system exhibits a considerable net magnetization.

### E_overflip.jpg Analysis

Unexpected observations from the graph include:

- Undoped system reaches equilibrium faster than lightly doped ones.
- Equilibrium time for 0.5% ~ 3% doping is similar, while 10% ~ 40% doping reaches equilibrium faster.

### DE_vs_C_log.jpg Analysis

The graph shows a roughly linear relationship between 0.5% ~ 40% doping, with the relationship equation:

\[ \Delta E \propto C^{0.6} \]

This suggests that as the concentration of disorder increases, the interaction block per site decreases.

## Validation

Validation is provided through:

- Examination of undoped Ising model behavior at low and high temperatures.
- Animation detailing flip dynamics through various concentrations of disorder.

## Future Exploration

Due to time limitations, certain questions remain unexplored. Future explorations could include studying the effects of different disorders, survival time length of magnetization, and the relationship between alignment time, concentration, and other factors.

## Note on Functionality

Every non-commented function should be operational, but they haven't been extensively tested on other machines. If issues arise, try deleting `.o` and `.x` files and recompile if necessary. Additionally, further testing on various computers is recommended.

Feel free to explore and expand upon this work to address the remaining questions and potential improvements.
