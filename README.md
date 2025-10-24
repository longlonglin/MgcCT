

# Environment Setup

	The proposed algorithms are implemented in Python 3.10.12. All experiments and performance evaluations reported in this paper are conducted on a high-memory Linux server equipped with a single Intel Xeon Silver 4210 processor (2.20GHz) and 1 TB of RAM.

	For development and testing purposes, our codebase maintains cross-platform compatibility. It has been verified to run on consumer-grade hardware, including an Apple M1 system with 8GB of RAM running macOS Monterey 12.3, as well as an Intel Core i7-10700 system with 16GB of RAM running Windows 10.

	However, for processing large-scale datasets and achieving optimal performance, we strongly recommend utilizing a server-grade environment. The primary server configuration used in our experiments provides the necessary computational resources—particularly the ample memory capacity—to ensure efficient execution and stability.





 # Usage

python mgc.py cm.txt

The  running results are as follows

	cm.txt is loading...
	loading_graph_time(s)0.15706634521484375
	Tspan:193.7057986111111
	number of static edges: 20292
	number of temporal edges: 59794
	number of nodes: 1893
	################./motifs/M1.txt###########
	delta3600.0
	################AGA###########
	TMC_quality_AGA  0.008264462809917356
	time of AGA  1.8325653076171875
	################HLA###########
	time of HLA0.0026219844818115234
	TMC_quality_HLA0.03511734164023254
	################./motifs/M2.txt###########
	delta3600.0
	################AGA###########
	TMC_quality_AGA  0.0014144271570014145
	time of AGA  1.1148264408111572
	################HLA###########
	time of HLA0.0006194591522216796
	TMC_quality_HLA0.0017775956259844304
	delta7200.0
	################AGA###########
	TMC_quality_AGA  0.005714285714285714
	time of AGA  1.2868828773498535
	################HLA###########
	time of HLA0.0016592979431152345
	TMC_quality_HLA0.021860131101142592
	delta10800.0
	################AGA###########
	TMC_quality_AGA  0.0047169811320754715
	time of AGA  1.4242990016937256
	################HLA###########
	time of HLA0.0007436752319335937
	TMC_quality_HLA0.020650517045459886
	delta14400.0
	################AGA###########
	TMC_quality_AGA  0.0045871559633027525
	time of AGA  1.5636961460113525
	################HLA###########
	time of HLA0.0014356136322021484
	TMC_quality_HLA0.0599399093436737
	delta18000.0
	################AGA###########
	TMC_quality_AGA  0.0007220216606498195
	time of AGA  1.6116528511047363
	################HLA###########
	time of HLA0.0016808032989501954
	TMC_quality_HLA0.0555195777064301
	################./motifs/M3.txt###########
	delta3600.0
	################AGA###########
	TMC_quality_AGA  0.45492813680235
	time of AGA  1.7497291564941406
	################HLA###########
	time of HLA0.0019198417663574218
	TMC_quality_HLA0.6923997830575496
	################./motifs/M4.txt###########
	delta3600.0
	################AGA###########
	TMC_quality_AGA  0.00851581508515815
	time of AGA  1.8004209995269775
	################HLA###########
	time of HLA0.0006944656372070313
	TMC_quality_HLA0.06795732450791721
	################./motifs/M5.txt###########
	delta3600.0
	################AGA###########
	TMC_quality_AGA  0.43507972665148065
	time of AGA  5.139225482940674
	################HLA###########
	time of HLA0.0010783672332763672
	TMC_quality_HLA0.6364895402284989
	################./motifs/M6.txt###########
	delta3600.0
	################AGA###########
	TMC_quality_AGA  0.6153846153846154
	time of AGA  4.906836271286011
	################HLA###########
	time of HLA0.0009699821472167969
	TMC_quality_HLA0.8835071415679774

 
 

