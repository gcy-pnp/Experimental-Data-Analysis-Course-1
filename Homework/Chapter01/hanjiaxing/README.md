# HJX's homework 01 #

* 1st Step: Simulation
  * root sim.cpp，生成 sim.root

* 2nd Step: X Position Calibration
  * root pos1.cpp，微分法确定 td-tu 均匀分布的边界，生成 pos1.root
  * root pos2.cpp，刻度中子或gamma打在探测器的位置，生成 pos2.root

* 3rd step: TOF and Energy Calibration
  * root ene1.cpp, 确定gamma打在探测器中间的TOF，生成 ene1.root
  * root ene2.cpp，中子TOF的距离修正，得到中子能谱，生成 ene2.root
