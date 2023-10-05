# Projecting "Computational Modeling of a Biological Network Implementing the Rescorla-Wagner Learning Rule: A Simulation Study"
The article in the reprository describe all steps, you could follow code base on that. Different codes for each section is available as .m file in the the folder "MATLAB"

This repository contains code and visualizations related to the simulation of neural pathways involving Projecting Neurons (PNs), Mushroom Body Input Neurons (MBIN), Kenyon Cells (KCs), and Mushroom Body Output Neurons (MBON) within a computational model. The code demonstrates how inputs from different sensory stimuli are processed through these neural components.

## Prerequisites
Ensure you have the following libraries installed before running the code:
- `matplotlib`
- `numpy`

## Getting Started
To visualize the neural pathways, execute the provided Python code in your local environment or a Jupyter notebook. Make sure to have the required libraries installed.

```python
import matplotlib.pyplot as plt
import torch
import numpy as np
import matplotlib.gridspec as gridspec
#globals().clear()  # Clears all variables in the global scope

# ... [The code you provided for data generation and processing]

# Adjust spacing between subplots
plt.tight_layout()

# Show the plot
plt.show()
```

## Understanding the Code
1. **Inputs to the System:** The code initializes different sensory inputs (CS_1, CS_2, US_1, US_2) and combines them to form the input signals (CS1, CS2, US11, US12, US21, US22).

2. **Kenyon Cells (KCs):** The code models the activity of Kenyon Cells using a combination of inputs and an activation function. The KC activities are visualized in separate subplots.

3. **Mushroom Body Output Neurons (MBONs):** The code simulates the Mushroom Body Output Neurons by considering inhibitory and excitatory projections. The synaptic plasticity (STDP) is applied to update synaptic weights based on the neural activities.

## Visualization
The code generates visualizations for:
- **Inputs to KCs:** Demonstrated through subplots displaying CS (Odour) and US activities to first and second KCs.
- **KC Activities:** Visualized for the first and second KCs.
- **MBON Activities:** Displayed for Approach (M1) and Avoid (M2) neurons.

## Additional Notes
- Ensure that the input data and parameters are correctly set for meaningful simulation results.
- Experiment with different input patterns and parameters to observe the behavior of the neural network.

Feel free to modify and extend the code as needed for your experiments. For any questions or issues, please refer to the provided documentation or contact the repository owner. Happy coding!
