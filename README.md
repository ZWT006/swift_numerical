# kinodynamicpath

#### 介绍
HUST AIA Coralab Aggressive Navigation MATLAB code
应用于Quadruped Robot的 Fast Realtime Navigation under Unknown Environment 

#### 代码结构
```
|-- LazyKinoPRM.m
|-- LazyPRM
|   |-- AngleCost.m
|   |-- AngleDelta.m
|   |-- CloseList.m
|   |-- CollisionCheck.m
|   |-- DistanceCost.m
|   |-- ExpandNode.m
|   |-- ObstacleCheck.m
|   |-- OpenList.m
|   |-- Point.m
|   |-- TrajectoryCheck.m
|   `-- insertGoalNode.m
|-- LocalOBVP
|   |-- AngleDelta.m
|   |-- BVP.m
|   |-- OBVP.m
|   |-- SampleOBVP.m
|   |-- formula_compute.mlx
|   `-- obvp_compute.mlx
|-- PolyOpt
|   |-- CostFunc.m
|   |-- PolyTraj.m
|   |-- calculate.mlx
|   |-- dynamicCost.m
|   |-- getAbeqMatrix.m
|   |-- getCoeffCons.m
|   |-- matlabopt.m
|   |-- memorandum.txt
|   |-- nloptTO.m
|   |-- nonConstrain.m
|   |-- obstacleCost.m
|   |-- optCoeffs.xlsx
|   |-- ovalCost.m
|   |-- sdfMap.m
|   |-- sdfMapGen.m
|   |-- smoothCost.m
|   |-- test.m
|   `-- timeCost.m
|-- README.en.md
|-- README.md
|-- TrajGener
|   |-- BerneteinCoeff.m
|   |-- BezierCurve.m
|   |-- CartesianPolynomials.m
|   |-- ChainedForm.m
|   |-- MinimumPolySolver.m
|   |-- SegPolyTraj.m
|   |-- VelocityVector.m
|   |-- calcu.mlx
|   |-- calcutest.m
|   |-- getAMbeq.m
|   |-- getAMbieq.m
|   |-- getAbeq.m
|   |-- getBAbeq.m
|   |-- getBAbieq.m
|   |-- getCoeffCons.m
|   |-- getCt.m
|   |-- getM.m
|   |-- getMQM.m
|   |-- getQ.m
|   |-- hw1_1.m
|   |-- hw1_2.m
|   |-- nearestSPD.m
|   |-- picture
|   `-- plotaxis.m
|-- TrajOpt
|   |-- BerneteinCoeff.m
|   |-- BezierCurve.m
|   |-- SegPolyTraj.m
|   |-- TrajectoryOpt.m
|   |-- getAbeq.m
|   |-- getBAbeq.m
|   |-- getBAbieq.m
|   |-- getCoeffCons.m
|   |-- getM.m
|   |-- getMQM.m
|   |-- getQ.m
|   |-- maptest.m
|   |-- myconstraint.m
|   |-- myfunc.m
|   |-- nearestSPD.m
|   |-- nloptTO.m
|   |-- path.mat
|   |-- sdf_map.png
|   |-- symscalcu.mlx
|   `-- test.m
|-- commond.txt
|-- map
|   |-- map0.png
|   |-- map09.png
|   |-- map1.png
|   |-- map2.png
|   |-- map3.png
|   |-- map4.png
|   |-- map5.png
|   |-- map6.png
|   |-- map7.png
|   |-- map8.png
|   `-- map9.png
|-- opt_fig.m
|-- path.mat
|-- pose_map.mat
|-- test.m
|-- test_a_obvp.m
|-- testblog.txt
`-- traj_test.m
```


#### 代码介绍

1.  `LazyPRM`   :LazyPRM 算法实现,主要包括路径搜索的几个关键函数以及`OpenList.m`,`CloseList.m`,`Point.m`
2.  `LocalOBVP` :LocalOBVP 路径搜索过程中使用 Optimal Boundary Value Problem 来连接Node,主要使用 `SampleOBVP.m`
3.  `map`       :保存的障碍地图
4.  `PolyOpt`   :Nonlinear Trajectory Optimization,主要文件为`matlatopt.m`,进行非线性迭代优化
5.  `TrajGener` :分段曲线优化,主要文件为`SegPolyTraj.m`,用于**多项式**以及**贝塞尔**曲线优化
6.  `TrajOpt`   :内容与`TrajGener`类似,不用care

#### 使用说明

1.  初始路径搜索：`LazyKinoPRM.m`进行路径初始搜索,相关参数见具体文件
2.  轨迹分段优化：`PolyOpt/matlabopt.m`进行路径的分段非线性优化,相关参数见具体文件

#### 参与贡献

1.  Fork 本仓库
2.  新建 Feat_xxx 分支
3.  提交代码
4.  新建 Pull Request


#### 特技

1.  使用 Readme\_XXX.md 来支持不同的语言，例如 Readme\_en.md, Readme\_zh.md
2.  Gitee 官方博客 [blog.gitee.com](https://blog.gitee.com)
3.  你可以 [https://gitee.com/explore](https://gitee.com/explore) 这个地址来了解 Gitee 上的优秀开源项目
4.  [GVP](https://gitee.com/gvp) 全称是 Gitee 最有价值开源项目，是综合评定出的优秀开源项目
5.  Gitee 官方提供的使用手册 [https://gitee.com/help](https://gitee.com/help)
6.  Gitee 封面人物是一档用来展示 Gitee 会员风采的栏目 [https://gitee.com/gitee-stars/](https://gitee.com/gitee-stars/)
