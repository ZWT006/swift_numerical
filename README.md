# swift_numerical
**Huazhong University of Science and Technology**, **School of Artificial Intelligence and Automation**, **Co**ntrol **R**obotics and **A**utomation Laboratory (CoraLab) *Hierarchical Planner for Agile Quadruped Navigation with Awareness of Motion Anisotropy* for Numerical Calculation MATLAB code
#### Description
1. **LazyPRM**: The implementation of the LazyPRM algorithm, which mainly includes several key functions for path search as `OpenList.m`, `CloseList.m`, and `Point.m`.
2. **LocalOBVP**: LocalOBVP uses Optimal Boundary Value Problem to connect Nodes in the path search process, mainly using `SampleOBVP.m`.
3. **map**: store obstacle map.
4. **PolyOpt**: Nonlinear Trajectory Optimization, mainly implemented in `matlatopt.m` using nonlinear iterative optimization.
5. **TrajGener**: Piecewise curve optimization, mainly implemented in `SegPolyTraj.m` for polynomial and Bezier curve optimization.


#### Instructions
First, use `LazyKinoPRM.m` for initial path search. Then, use `matlabopt.m` for trajectory optimization to obtain a continuous and smooth spatio-temporal trajectory.
**Kinodynamic Trajectory Generation**
`mapaddress`: Path of the obstacle map.
`Point_i/Point_f`: Initial point / target point represented as $[x,y,θ]$, indicating horizontal displacement and orientation angle.
**Nonlinear Trajectory Optimization**
`seq_start/seq_end`: Starting and ending segment indices for optimization.

#### Postscript
1. The MATLAB code here is the theoretical verification part of our project. This algorithm is applied in practical navigation planning.
2. If you have any questions or suggestions, please feel free to contact `zwt190315@163.com` or `wentaozhang@hust.edu.cn`. Thank you very much.
