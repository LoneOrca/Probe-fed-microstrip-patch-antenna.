# Project: 16-GHz Probe-Fed Microstrip Patch Antenna

## Project Overview
This project involves the design, calculation, and simulation of a **16-GHz microstrip patch antenna** using **Ansys HFSS**. The design is implemented on a **Duroid 5880 substrate** with a probe-fed configuration.

---

## Technical Specifications

* **Operating Frequency ($f_0$):** 16 GHz 
* **Substrate:** Rogers Duroid 5880 ($\epsilon_r = 2.2$) 
* **Substrate Thickness ($h$):** 20 mils (0.508 mm) 
* **Feeding Method:** Probe-fed using a 0.085" semi-rigid coax cable 
* **Coax Dielectric:** Teflon ($\epsilon_r = 2.1$) 
---

## Design & Analysis Highlights

The project follows a rigorous workflow from theoretical calculations to full-wave EM simulation:

### 1. Mathematical Design

Theoretical parameters were calculated using MATLAB to determine the physical dimensions of the antenna.

* **Patch Length ($L_p$):** 5.78369 mm 
* **Patch Width ($W$):** 8.67554 mm 
* **Feed Position ($x_f$):** Optimized for 50 $\Omega$ impedance matching.

### 2. HFSS Simulation & Tuning

The antenna was modeled and tuned in **Ansys HFSS** to achieve high performance.

* **Input Return Loss:** Optimized to be better than 25 dB ($|S_{11}| \le -25$ dB) at 16 GHz.
* **Gain:** Achieved a simulated gain of **8.3886 dB**.
* **Bandwidth:** Demonstrated a 10 dB return loss bandwidth of approximately **3.88%** (625 MHz).

### 3. Radiation Characteristics

* **E-Plane Beamwidth ($BW_{EP}$):** 61.0790° 
* **H-Plane Beamwidth ($BW_{HP}$):** 79.0974° 

---

## Repository Contents

* **Simulation:** HFSS project files and design properties.

* **Calculations:** MATLAB scripts for design parameter verification.

* **Results:** Detailed tables comparing calculated vs. simulated data.

* **Plots:** $S_{11}$ frequency response and 2D/3D radiation patterns.

---

## Tools Used

* **Ansys HFSS 2023 R2 (Student Edition):** EM Simulation.

* **MATLAB:** Numerical calculations.
---
Developed as part of the Wave Transmission and Reception (EE 557/457) curriculum at Western New England University.
