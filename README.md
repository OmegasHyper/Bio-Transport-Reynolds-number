# Reynolds Number Flow Analysis (Biotransport SBEG201)

<img width="1280" height="853" alt="Image" src="https://github.com/user-attachments/assets/2aaf93bc-d708-48d0-b5ab-966704b41cc6" />

This repository contains the computational analysis of the **Reynolds Number ($\mathbf{Re}$)** and its role in classifying fluid flow into Laminar, Transitional, and Turbulent regimes. The project uses advanced software integration to fulfill analytical and visualization requirements.

***

## 🎯 Project Highlights & Core Achievements

1.  **Analytical Modeling:** Calculated $\mathbf{Re}$ for a comprehensive matrix of configurations, including critical biofluids (Blood, CSF), demonstrating the transition from viscous to inertial dominance.
2.  **Advanced Architecture:** Implemented a **hybrid Node.js/MATLAB system** for the interactive calculator, separating the front-end user experience from the core analytical engine.
3.  **CFD Simulation (Bonus):** Applied **OpenFOAM** to model high-$\mathbf{Re}$ turbulent flow (Lid-Driven Cavity), visualizing complex momentum transfer patterns.

***

## 🛠️ Technology Stack & Implementation

| Component | Tool / Language | Key Function |
| :---: | :---: | :--- |
| **Core Analysis & Visualization** | **MATLAB R2025a** | $\mathbf{Re}$ Calculation, Plotting of theoretical profiles ($\mathbf{f}$ vs. $\mathbf{Re}$, Velocity Profiles). |
| **Interactive Application** | **Node.js + MATLAB** | Server bridge for seamless, cross-platform communication between the GUI and the analytical engine. |
| **Advanced Simulation** | **OpenFOAM (CFD)** | Modeling of fluid dynamics under high inertial conditions. |

***

## ✨ Key Visual Results

### 1. Hybrid Calculator Verification
The application provides instant, color-coded regime classification.
<img width="1044" height="726" alt="Image" src="https://github.com/user-attachments/assets/90b903d7-2364-459b-950c-e88262847c4e" />
<img width="1065" height="729" alt="Image" src="https://github.com/user-attachments/assets/caa5a726-7278-44c6-9c66-e5cb10c7901d" />

### 2. Profile Transformation
Visualization shows the fundamental shift from the **Parabolic** profile (Laminar, $\mathbf{V_{\text{max}} = 2 V_{\text{avg}}}$) to the **Blunt** turbulent profile.
<img width="707" height="725" alt="Image" src="https://github.com/user-attachments/assets/4836f3e8-92e6-4738-967c-608fdbf2550c" />

### 3. Frictional Resistance Analysis
The logarithmic plot confirms the energy loss shift: steep decline in the laminar regime ($\mathbf{f \propto 1/Re}$) versus a shallow decline in the turbulent regime ($\mathbf{f \propto 1/Re^{0.25}}$).


<img width="725" height="731" alt="Image" src="https://github.com/user-attachments/assets/59c681d0-9c82-4b85-93f7-34f299f657f0" />


### 4. CFD Flow Visualization 
The OpenFOAM result illustrates the stable, dominant primary vortex formed by intense momentum transfer under high **Re** conditions.

<img width="1132" height="794" alt="Image" src="https://github.com/user-attachments/assets/6acdffa5-e56b-429b-8182-ba88706a7d1a" />
<img width="1132" height="795" alt="Image" src="https://github.com/user-attachments/assets/00610b43-f987-4fc6-b6a6-2fd66a731e9d" />

## 📁 Repository Structure

* `Assignment_Report.pdf`: The complete, detailed academic report.
* `MATLAB_Scripts/`: Contains all analytical `.m` files (Parts 1-5).
* `GUI_NodeJS/`: Files for the interactive application server and front-end.
* `CFD_OpenFOAM/`: Simulation case files and final visualization images.

---
**Developed by:** [Team 7]
- @OmegasHyper
- @Amatalrahman
- @Alaa-Essam5
- @engy27005
- @engyelsarta

**Course:** SBEG201 Biotransport
