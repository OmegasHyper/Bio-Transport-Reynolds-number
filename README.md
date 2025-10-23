# 🔬 Reynolds Number Flow Analysis (Biotransport SBEG201)

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

> **[Insert Image of GUI showing Green (Laminar) and Red (Turbulent) results here]**

### 2. Profile Transformation
Visualization shows the fundamental shift from the **Parabolic** profile (Laminar, $\mathbf{V_{\text{max}} = 2 V_{\text{avg}}}$) to the **Blunt** turbulent profile.

> **[Insert Image of Velocity Profiles vs Reynolds Number (Laminar vs. Turbulent) here]**

### 3. Frictional Resistance Analysis
The logarithmic plot confirms the energy loss shift: steep decline in the laminar regime ($\mathbf{f \propto 1/Re}$) versus a shallow decline in the turbulent regime ($\mathbf{f \propto 1/Re^{0.25}}$).

> **[Insert Image of Darcy Friction Factor vs Reynolds Number Plot here]**

### 4. CFD Flow Visualization (Bonus)
The OpenFOAM result illustrates the stable, dominant primary vortex formed by intense momentum transfer under high **Re** conditions.

> **[Insert Image of OpenFOAM/CFD Lid-Driven Cavity Flow result here]**

***

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
