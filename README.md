

# Environment Setup

Our algorithms are implemented in Python 3.7.3  and  all experiments are conducted on a Linux server with one Intel(R) Xeon(R) Silver 4210 2.20GHz CPU and 1 TB RAM.

Ps: Our code can also run on a desktop with Apple M1 and 8GB RAM running macOS Monterey 12.3  (and Inter(R) Core(TM) i7-10700@2.90GHZ and 16 GB RAM running Windows 10). But we recommend you run it on a server because the server has enough memory to handle large datasets and run faster.




 # Usage

python ECRC.py Cora.txt

The  running results are as follows

	Cora.txt is loading...
	number of left_nodes2708
	number of right_nodes1432
	average_left_degree18.174298375184637
	average_right_degree34.36871508379888
	number of edges49216.0
	p,q: 2  2
	################ECRC###########
	quality_ECRC  0.13350874452841222
	time of ECRC  40.219770193099976
	################ECRC_E###########
	quality_ECRC_E  0.01697685939889084
	time of ECRC_E  2.20005202293396
	p,q: 2  3
	################ECRC###########
	quality_ECRC  0.04286252203668915
	time of ECRC  40.4626841545105
	################ECRC_E###########
	quality_ECRC_E  0.0004645877621797628
	time of ECRC_E  2.203089952468872
	p,q: 3  2
	################ECRC###########
	quality_ECRC  0.0
	time of ECRC  31.39865517616272
	################ECRC_E###########
	quality_ECRC_E  0.004221901818285617
	time of ECRC_E  2.1270155906677246
	p,q: 2  4
	################ECRC###########
	quality_ECRC  0.017497699365793568
	time of ECRC  33.816262006759644
	################ECRC_E###########
	quality_ECRC_E  1.7294161220552104e-05
	time of ECRC_E  2.221113443374634
	p,q: 3  3
	################ECRC###########
	quality_ECRC  0.09766514471431044
	time of ECRC  31.932520389556885
	################ECRC_E###########
	quality_ECRC_E  1.585627553892494e-05
	time of ECRC_E  2.0966732501983643
	p,q: 4  2
	################ECRC###########
	quality_ECRC  0.24280150603998957
	time of ECRC  24.059733152389526
	################ECRC_E###########
	quality_ECRC_E  0.0013038401534577252
	time of ECRC_E  2.122148275375366

 
Our model has only two parameters, ``p`` and ``q``, which users can easily modify by changing line 519 of ECRC.py.
 

