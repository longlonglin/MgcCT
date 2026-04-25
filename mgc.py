import time
import sys
from collections import defaultdict
from heapq import *
from random import choice
from binary_heap import *
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.csgraph import laplacian
import itertools
from TM import *

class Graph:

	def __init__(self, dataset):
		self.tadj_list = self.TemporalGraph(dataset)
		s_edge, t_edge, self.Tspan, self.node = 0, 0,0,set()
		self.G=[]
		self.G_detemporal=[]
		tmp=set()
		for u in self.tadj_list:
			s_edge += len(self.tadj_list[u])
			self.node.add(u)
			for v in self.tadj_list[u]:
				self.node.add(v)
				t_edge += len(self.tadj_list[u][v])
				for t in self.tadj_list[u][v]:
					tmp.add(t)
					self.G.append((u,v,t))
				self.G_detemporal.append((u,v,1))
		self.Tspan=(max(tmp)-min(tmp))/86400 #days
		print("Tspan:" + str(self.Tspan))
		print("number of static edges: " + str(s_edge))
		print("number of temporal edges: " + str(t_edge))
		print("number of nodes: " + str(len(self.node)))
		self.t_edge=t_edge
		

	# input graph is directed by default
	def TemporalGraph(self,dataset):
		tadj_list = {}
		tmp=defaultdict(set) #undirected static graph for identify the largest weakly connected component (cc) 
		print(dataset + " is loading...")
		starttime = time.time()
		with open(dataset, 'r') as file:
			lines = file.readlines()
			for line in lines:
				line = line.split()
				f_id, t_id, time_id = int(line[0]), int(line[1]), float(line[2])
				if f_id == t_id: #delete the self-loop
					continue
				if f_id in tadj_list:
					if t_id in tadj_list[f_id]:
						tadj_list[f_id][t_id].add(time_id)
					else:
						tadj_list[f_id][t_id] = {time_id}
				else:
					tadj_list[f_id] = {}
					tadj_list[f_id][t_id] = {time_id}
				tmp[f_id].add(t_id)
				tmp[t_id].add(f_id)
		# find the largest cc		
		visited = set()
		largest_cc = set()
		for node in tmp:
			if node not in visited:
				cc = set()
				visited.add(node)
				Q = [node]
				while Q:
					v = Q.pop()
					cc.add(v)
					for u in tmp[v]:
						if	u not in visited:
							visited.add(u)
							Q.append(u)
				if len(cc)>len(largest_cc):
					largest_cc=cc
		#find the real input graph			
		for u in set(tadj_list):
			if u not in largest_cc:
				del tadj_list[u]
				continue
			for v in set(tadj_list[u]):
				if v not in largest_cc:
					del tadj_list[u][v]
				if len(tadj_list[u])==0:
					del tadj_list[u]
		endtime = time.time()
		print("loading_graph_time(s)" + str(endtime - starttime))
		return tadj_list

	def weighted_for_temporal(self,M,delta):
		self.weighted_graph=defaultdict(lambda: defaultdict(int))
		self.edge_weighted_graph=defaultdict(lambda: defaultdict(int))
		self.M_vertex_number=0
		tmp=set()
		for u,v in M:
			tmp.add(u)
			tmp.add(v)
		self.M_vertex_number=len(tmp)
		starttime=time.time()
		counter = OptimizedTemporalSubgraphCounter(self.G, M, delta)
		edge_counts,self.matches = counter.count()
		for (u,v) in edge_counts:
			if (v,u) not in edge_counts:
				self.edge_weighted_graph[u][v]=edge_counts[(u,v)]
				self.edge_weighted_graph[v][u]=edge_counts[(u,v)]
			else:
				self.edge_weighted_graph[u][v]=edge_counts[(u,v)]+edge_counts[(v,u)]
				self.edge_weighted_graph[v][u]=edge_counts[(u,v)]+edge_counts[(v,u)]
		runtime_edge_weighted_graph=time.time()-starttime
		return runtime_edge_weighted_graph



	#build weighted graph based on temporal motif
	def build_graph(self):
		weighted_graph=self.edge_weighted_graph
		#find maximum connected componet from the  weighted_graph, which will be used for clustering
		visited = set()
		largest_cc = set()
		for node in weighted_graph:
			if node not in visited:
				cc = set()
				visited.add(node)
				Q = [node]
				while Q:
					v = Q.pop()
					cc.add(v)
					for u in weighted_graph[v]:
						if  u not in visited:
							visited.add(u)
							Q.append(u)
			if len(cc)>len(largest_cc):
				largest_cc=cc
		return weighted_graph, largest_cc


	def cut_vol(self,S):
		vol=0
		cut=0
		motif_vol=0
		for i, match in enumerate(self.matches):
			tmp=set()
			for j, (u, v, t) in enumerate(match):
				tmp.add(u)
				tmp.add(v)
			tmp_len=len(tmp.intersection(S))
			vol+=tmp_len
			motif_vol+=self.M_vertex_number
			if 1<=tmp_len<=self.M_vertex_number-1:
				cut+=1
		return cut,vol,motif_vol
		
		


	
	def AGA(self): 
		starttime=time.time()
		motif_weighted_graph, motif_max_cc = self.build_graph()
		old_new_id = {}
		new_old_id = {}
		index = 0
		for u in sorted(list(motif_max_cc)):
			old_new_id[u]=index
			new_old_id[index]=u
			index+=1
		weighted_sum=0
		n=len(motif_max_cc)
		row=[]
		col=[]
		new_motif_weighted_graph,new_degree_motif={},{}
		for u in motif_max_cc:
			new_degree_motif[old_new_id[u]]=0
			new_motif_weighted_graph[old_new_id[u]]={}
			for v in motif_weighted_graph[u]:
				row.append(old_new_id[u])
				col.append(old_new_id[v])
				new_motif_weighted_graph[old_new_id[u]][old_new_id[v]]=motif_weighted_graph[u][v]
				new_degree_motif[old_new_id[u]]+=motif_weighted_graph[u][v]
				weighted_sum +=motif_weighted_graph[u][v]
		A = sp.csr_matrix((np.ones(len(row)), (row, col)), shape=(n, n))
		norm_L = laplacian(A,normed=True)
		emb_eig_val, p = spla.eigsh(norm_L, k=2, which='SM')
		pi = np.real(p[:, 1])
		pi = np.argsort(pi)
		S = set()
		volS = 0
		cutS = 0
		best_condu, best_index, count = float("inf"), 0, 0
		best_set = set()
		for x in pi:
			u = x
			S.add(u)
			count += 1
			for node in new_motif_weighted_graph[u]:
				if node in S:
					cutS -= 2 * new_motif_weighted_graph[u][node]
			cutS = cutS + new_degree_motif[u]
			volS = volS + new_degree_motif[u]
			if min(volS,  weighted_sum - volS) != 0 and cutS / min(volS,  weighted_sum - volS) < best_condu:
				best_condu = cutS / min(volS,  weighted_sum - volS)
				best_index = count
		for x in range(best_index):
			best_set.add(pi[x])
		if len(best_set) > len(new_motif_weighted_graph) / 2:
			best_set = set(new_motif_weighted_graph) - set(best_set)
		endtime=time.time()
		runtime=endtime-starttime
		#over and next compute higher-order conductance
		S=set()
		for v in best_set:
			S.add(new_old_id[v])	
		cut,vol,motif_vol=self.cut_vol(S)
		if min(vol,motif_vol-vol)==0:
			ans=1
		else:
			ans=cut/min(vol,motif_vol-vol)
		return ans, S, runtime,motif_max_cc


	def HLA(self,alpha,seed,epsilon):
		starttime=time.time()
		pi, r = {}, {}
		r[seed] = 1
		q = [seed]
		D=set()
		dm={}
		wef={}
		while q:
			u = q.pop()
			if u not in D:
				D.add(u)
				dm[u]=0
				wef[u]={}
				for v in self.edge_weighted_graph[u]:
					if self.edge_weighted_graph[u][v]!=0:
						wef[u][v]=self.edge_weighted_graph[u][v]
						dm[u]+=wef[u][v]
			if r[u]>=epsilon*dm[u]:
				for v in wef[u]:
					if v not in r:
						r[v] = 0
					r[v] = r[v] + (1 - alpha) * r[u] * wef[u][v] / dm[u]
					q.append(v)
				if u not in pi:
					pi[u] = 0
				pi[u] = pi[u] + alpha * r[u]
				r[u] = 0
		Sort_pi={}
		for u in pi:
			if pi[u]!=0 and dm[u]!=0:
				Sort_pi[u] = pi[u] / dm[u]
		pi = sorted(Sort_pi.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
		S = set()
		volS = 0
		cutS = 0
		best_condu, best_index, count = float("inf"), 0, 0
		best_set = set()
		for x in pi:
			u = x[0]
			S.add(u)
			count += 1
			for node in wef[u]:
				if node in S:
					cutS -= 2 * wef[u][node]
			cutS = cutS + dm[u]
			volS = volS + dm[u]
			if volS!=0 and cutS /volS < best_condu:
				best_condu = cutS /volS
				best_index = count
		for x in range(best_index):
			best_set.add(pi[x][0])
		if len(best_set) > len(self.tadj_list) / 2:
			best_set = set(self.tadj_list) - set(best_set)
		endtime=time.time()
		runtime=endtime-starttime
		
		#over and next compute higher-order conductance
		S=set()
		for v in best_set:
			S.add(v)
		cut,vol,motif_vol=self.cut_vol(S)
		if min(vol,motif_vol-vol)==0:
			ans=1
		else:
			ans=cut/min(vol,motif_vol-vol)
		return ans, S, runtime

		
if __name__ == "__main__":
	dataset = sys.argv[1]
	G = Graph(dataset)	
	t=86400/24

	motif_configs = {
		1: [t],
		2: [t,2*t,3*t,4*t,5*t],
		3: [t],
		4: [t],
		5: [t],
		6: [t]
	}

	for i, deltas in motif_configs.items():
		mpath = f"./motifs/M{i}.txt"
		print(f"################{mpath}###########")
		
		M = []
		with open(mpath, 'r') as file:
			lines = file.readlines()
			for line in lines:
				line = line.split()
				M.append((str(line[0]), str(line[1])))
		for delta in deltas:
			print(f"delta{delta}")
			t = G.weighted_for_temporal(M, delta)
			print(f"################AGA###########")
			quality, S, runtime,motif_max_cc=G.AGA()
			print("TMC_quality_AGA  "+str(quality))
			print("time of AGA  "+str(t+runtime))

			print(f"################HLA###########")
			alpha=0.2
			epsilon=1/G.t_edge
			avg_time = 0
			avg_quality = 0
			avg_len=0
			seed_number=50
			for i in range( seed_number):
				seed=int(choice(list(motif_max_cc)))
				quality,S,runtime= G.HLA(alpha,seed,epsilon)
				avg_time += runtime
				avg_quality += quality
			print("time of HLA" + str(avg_time/ seed_number))
			print("TMC_quality_HLA" + str(avg_quality / seed_number))

		
