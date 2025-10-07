# Damping_Redistribution
Data and code of the article "Impact of Angle-Voltage Coupling on Small-Signal Stability of Power Systems: A Damping Perspective"

This source code repository accompanies the following paper:

> X. Peng, C. Fu, P. Yang, X. Ru, and F. Liu, "Impact of Angle-Voltage Coupling on Small-Signal Stability of Power Systems: A Damping Perspective," *IEEE Transactions on Circuits and Systems I: Regular Papers*, 2025.

The full paper and the source code can be found at: https://github.com/lingo01/Damping_Redistribution

For any questions or uses of the source codes, please feel free to contact the first author, Xiaoyu Peng (pengxy19@tsinghua.org.cn), and the corresponding author, Feng Liu (lfeng@mail.tsinghua.edu.cn).

**CITATION**: If you use this code in your work, whether directly or indirectly, please cite the above paper.

**LICENSE**: This work is licensed under the MIT License. See the [LICENSE](LICENSE) file in the repository for details.



# Introduction to the source code

All code is written in `MATLAB` and can be run directly in `MATLAB 2018a` or later versions.

The simulations are conducted on a modified IEEE 39-bus system with heterogeneous inverter-interface devices. 

* The network settings are borrowed from `MATPOWER` dataset, and the admittance matrix is stored as `.\data&figure\Network_netY_IEEE39Hete.mat` and employed in the simulations. 
* The device dynamic forms have been introduced in the paper, and their parameters are stored in `.\data&figure\ElementList_IEEE39Hete.xlsx` with a mapping relation defined in `settingElement.m`. Throughout the codes, the abbreviation Gen, CD and QD represent synchronous generators, conventional droop-controlled inverters and quadratic droop-controlled inverters, respectively.

The simulations are shown in Figs.8, 9, 10, and 11 in the paper. Next, we explain how to generate all the data and figures of this paper in detail.



## Fig.8 of this paper:

Run `main_Eig.m` to generate Fig.8 of this paper, where

* `main_EigLoad` generates Fig.8(a)-(d) corresponding to the eigenvalue changes under load scaling factors.
* `main_EigDevCtrl` generates Fig.8(e)-(h) corresponding to the eigenvalue changes under device voltage control strength scaling factors.
* `main_EigDevTime` generates Fig.8(i)-(l) corresponding to the eigenvalue changes under device voltage time constant scaling factors.



## Fig.9 of this paper:

Run `main_Det.m` to generate Fig.9 of this paper, where

* `main_DetLoad` generates Fig.9(a) showing the determinant ratio under different load scaling factors.
* `main_DetDevCtrl` generates Fig.9(b) showing the determinant ratio under different device voltage control strength scaling factors.
* `main_DetDevTime` generates Fig.9(c) showing the determinant ratio under different device voltage time constant scaling factors.
* `main_DetNetwork` generates Fig.9(d) showing the determinant ratio under different network connection strength scaling factors.
* `main_DetCtrlStrategy` generates Fig.9(e) showing the determinant ratio under different GFM-GFL penetration ratios.
* `main_DetWithEig` generates Fig.9(f)-(g) showing the determinant ratio with the relation to the real parts of the rightmost eigenvalues under different scenarios.



## Fig.10 of this paper:

It is noted that this figure is time-consuming to generate due to the low efficiency in solving the equilibrium of high-dimension and nonlinear power system models.

* `main_Dimension_uniform.m` generates the total damping and average damping under the situations of *uniform* distribution of integrated devices in Fig.10(a)-(b).
* `main_Dimension.m` generates the total damping and average damping under the situations of *nonuniform* distribution of integrated devices in Fig.10(a)-(b).
* `main_Dimension_Det_uniform.m` generates the determinant ratios under the situations of *uniform* distribution of integrated devices in Fig.10(c).
* `main_Dimension_Det.m` generates the determinant ratios under the situations of *nonuniform* distribution of integrated devices in Fig.10(c).
* The eigenvalue distributions displayed in Fig.10(d)-(e) are also generated in `main_Dimension_Det_uniform.m`.



## Fig.11 of this paper

Run `main_ControlEig.m` to generate Fig.11 of this paper.
