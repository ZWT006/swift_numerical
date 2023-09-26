# kinodynamicpath
**华中科技大学**，**人工智能与自动化学院**，**控制**,**机器人**与**自动化**实验室（CoraLab） *Hierarchical Planner for Agile Quadruped Navigation with Awareness of Motion Anisotropy*数值计算MATLAB代码

#### 代码介绍

1.  `LazyPRM`   :LazyPRM 算法实现,主要包括路径搜索的几个关键函数以及`OpenList.m`,`CloseList.m`,`Point.m`
2.  `LocalOBVP` :LocalOBVP 路径搜索过程中使用 Optimal Boundary Value Problem 来连接Node,主要使用 `SampleOBVP.m`
3.  `map`       :保存的障碍地图
4.  `PolyOpt`   :Nonlinear Trajectory Optimization,主要文件为`matlatopt.m`,进行非线性迭代优化
5.  `TrajGener` :分段曲线优化,主要文件为`SegPolyTraj.m`,用于**多项式**以及**贝塞尔**曲线优化
6.  `TrajOpt`   :内容与`TrajGener`类似,不用关心


#### 使用说明
首先使用`LazyKinoPRM.m`进行路径初始搜索，然后使用`matlabopt.m`进行轨迹优化，得到连续光滑的时空轨迹
初始路径搜索：`LazyKinoPRM.m`
这是一个参考使用LazyPRM*的轨迹生成方法，在空间中进行带有随机偏差的均匀采样，参数说明如下：
`mapaddress`: 障碍物地图路径
`Point_i/Point_f`: 起始点/目标点 $[x,y,\theta]$ 代表水平位移和方向角

轨迹分段优化：`matlabopt.m`
在多项式轨迹的基础上，进行轨迹的分段非线性优化，主要参数说明如下：
`seq_start/seq_end`: 优化的起始段序号

#### 备注
1. 这里的MATLAB代码是我们课题的理论验证部分，该算法应用在实际的导航规划中
2. 如果有疑问，欢迎联系`zwt190315@163.com`或者`wentaozhang@hust.edu.cn`提出建议，非常感谢

