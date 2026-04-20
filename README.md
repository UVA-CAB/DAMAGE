# DAMAGE
This repository provides an implementation of the DAMAGE (Diffuse Axonal Multi-Axis General Evaluation) metric, a physics-based model for estimating maximum brain strain (MPS) from rotational head kinematics.

The method is based on the paper:
“Development of a Second-Order System for Rapid Estimation of Maximum Brain Strain” 
(Gabler, L.F., Crandall, J.R. & Panzer, M.B. Ann Biomed Eng, 2018. https://doi.org/10.1007/s10439-018-02179-9)

Administrative contact: Ramakrishnan Iyer (zdf5ks@virginia.edu)

Research contact: Matthew B. Panzer (panzer@virginia.edu)

Repository Structure

├── LICENSE

├── damage.py                  # DAMAGE metric computation

│── example_data.xslx          # Example Head angular acceleration time history (Time (s), Angular Acceleration (rad/s/s))

└── README.md

Input

time = [t0, t1, ..., tn]           # units: s

r_ax = [...]                       # units: rad/s/s

r_ay = [...]                       # units: rad/s/s

r_az = [...]                       # units: rad/s/s

Ouput

DAMAGE
