# **Erbium-Doped Fiber Amplifier (EDFA) Matlab Model**

## **1. Project Overview**
This project involves the design and implementation of a **Matlab model** for **Erbium-Doped Fiber Amplifiers (EDFAs)** operating under static conditions. EDFAs are essential for loss compensation in fiber-optic communication systems, particularly in the communication spectrum based on silica fibers.

The project focuses on the development of a **comprehensive Matlab model** and generates various plots that provide insights into the performance of the EDFA.

## **2. Project Structure**
The project contains the following main components:
- **Matlab Simulation Files**: Contains the core `.m` files for simulating the behavior of EDFAs under different conditions.
- **Plots & Comments**: Various plots are generated from the simulation and included in the PDF report with comments.
- **Report**: A detailed PDF document that explains the Matlab model, plots, and analysis.

## **3. Deliverables**
The project deliverables include:
1. A **PDF report** containing:
   - Detailed Matlab model explanation.
   - The following **seven plots** extracted from simulation:
     - **a)** Absorption and emission cross-sections as extracted from slides using a plot digitizer.
     - **b)** \( P_p \) and \( P_{sout} \) vs position in the fiber (z) at a wavelength of **1550nm**.
     - **c)** **Gain vs fiber length** at different pump powers, including calculation of the **optimum fiber length** at **0.5W pump power**.
     - **d)** \( P_{sout} \) vs \( P_{sin} \) on a **log scale**, comparing saturation input power to theoretical values at the optimum fiber length.
     - **e)** **Gain vs \( P_{sin} \)** on a **log scale**, comparing saturation input power to theoretical values at the optimum fiber length.
     - **f)** **Gain vs wavelength (\( \lambda \))** at a certain pump power and input power.
     - **g)** \( P_{sout} \) with **ASE noise vs wavelength** on a log scale with input power at **1550nm**. The noise bandwidth is set to **0.1nm**.
   
2. A **zip file** containing:
   - The **Matlab simulation files (.m)** that generate the results and plots.

## **4. Key Project Specifications**
- **Wavelength of Operation**: **1550nm**.
- **Pump Power**: Simulation includes varying pump powers, with an emphasis on **0.5W** to calculate the optimum fiber length.
- **Gain Calculation**: Gain is calculated for different fiber lengths, pump powers, and input powers.
- **Noise Calculations**: ASE noise is calculated across a wavelength range with a **0.1nm increment** in the wavelength vector.
