import bisect
import sys
from collections import defaultdict
import time
import heapq
from itertools import islice

class OptimizedTemporalSubgraphCounter:
	def __init__(self, G_edges, M_edges, delta):


		self.G_edges = sorted(G_edges, key=lambda x: x[2])
		self.M_edges = M_edges
		self.delta = delta
		self.n = len(self.G_edges)
		

		self._build_indices()
		

		self._initialize_state()
		

		self.edge_counter = defaultdict(int)
		

		self.timestamps_array = self.timestamps  #

	def _build_indices(self):
		self.out_edges = defaultdict(list)

		self.in_edges = defaultdict(list)

		self.pair_edges = defaultdict(list)
		self.timestamps = []  
		

		for idx, (u, v, t) in enumerate(self.G_edges):
			self.timestamps.append(t)
			

			self.out_edges[u].append((t, v, idx))
			

			self.in_edges[v].append((t, u, idx))
			

			self.pair_edges[(u, v)].append((t, idx))
		

	def _initialize_state(self):

		self.map_GM = {} 
		self.map_MG = {}  
		

		G_nodes = set()
		for u, v, _ in self.G_edges:
			G_nodes.add(u)
			G_nodes.add(v)
		for node in G_nodes:
			self.map_GM[node] = -1

		M_nodes = set()
		for u, v in self.M_edges:
			M_nodes.add(u)
			M_nodes.add(v)
		for node in M_nodes:
			self.map_MG[node] = -1
		

		self.edge_count = defaultdict(int)
		

		self.e_stack = []  
		self.eG = 0	  
		self.eM = 0	  
		self.t0 = None	
	
	def _find_next_match(self):

		uM, vM = self.M_edges[self.eM]
		uG_mapped = self.map_MG.get(uM, -1)
		vG_mapped = self.map_MG.get(vM, -1)
		

		map_GM = self.map_GM
		G_edges = self.G_edges
		

		if uG_mapped != -1 and vG_mapped != -1:

			candidates = self._get_pair_candidates(uG_mapped, vG_mapped)
		elif uG_mapped != -1:

			candidates = self._get_out_candidates(uG_mapped)
		elif vG_mapped != -1:

			candidates = self._get_in_candidates(vG_mapped)
		else:

			candidates = self._get_global_candidates()
		

		for cand_idx in candidates:
			uG, vG, t = G_edges[cand_idx]
			
			u_ok = (uG_mapped != -1 and uG == uG_mapped) or \
				   (uG_mapped == -1 and map_GM.get(uG, -1) == -1)
			

			v_ok = (vG_mapped != -1 and vG == vG_mapped) or \
				   (vG_mapped == -1 and map_GM.get(vG, -1) == -1)
			
			if u_ok and v_ok:
				return cand_idx
		
		return self.n  

	def _get_pair_candidates(self, uG, vG):

		key = (uG, vG)
		if key not in self.pair_edges:
			return
		edges_tmp = self.pair_edges[key]
		t0 = self.t0
		eG = self.eG
		

		if t0 is not None:

			pos = bisect.bisect_right(edges_tmp, (t0, float('inf')))
			if pos == 0:
				return

			edges = edges_tmp[:pos]
		else:
			edges = edges_tmp

		lo, hi = 0, len(edges)
		while lo < hi:
			mid = (lo + hi) // 2
			if edges[mid][1] < eG:
				lo = mid + 1
			else:
				hi = mid
		

		for i in range(lo, len(edges)):
			yield edges[i][1]

	def _get_out_candidates(self, uG):

		if uG not in self.out_edges:
			return
		
		edges_tmp = self.out_edges[uG]
		t0 = self.t0
		eG = self.eG
		

		if t0 is not None:

			dummy = (t0, float('inf'))
			pos = bisect.bisect_right(edges_tmp, dummy)
			if pos == 0:
				return
			edges = edges_tmp[:pos]
		else:
			edges = edges_tmp
		

		lo, hi = 0, len(edges)
		while lo < hi:
			mid = (lo + hi) // 2
			if edges[mid][2] < eG:
				lo = mid + 1
			else:
				hi = mid
		

		for i in range(lo, len(edges)):
			yield edges[i][2]
	
	def _get_in_candidates(self, vG):

		if vG not in self.in_edges:
			return
		
		edges_tmp = self.in_edges[vG]
		t0 = self.t0
		eG = self.eG
		

		if t0 is not None:

			dummy = (t0, float('inf'))
			pos = bisect.bisect_right(edges_tmp, dummy)
			if pos == 0:
				return
			edges = edges_tmp[:pos]
		else:
			edges = edges_tmp
		

		lo, hi = 0, len(edges)
		while lo < hi:
			mid = (lo + hi) // 2
			if edges[mid][2] < eG:
				lo = mid + 1
			else:
				hi = mid
		

		for i in range(lo, len(edges)):
			yield edges[i][2]
	
	def _get_global_candidates(self):

		eG = self.eG
		n = self.n
		t0 = self.t0		

		if t0 is not None:

			end_pos = bisect.bisect_right(self.timestamps_array, t0, lo=eG)
			for idx in range(eG, end_pos):
				yield idx
		else:
			for idx in range(eG, n):
				yield idx
	
	def _update_mappings(self, edge_idx):

		uG, vG, _ = self.G_edges[edge_idx]
		uM, vM = self.M_edges[self.eM]
		

		self.map_GM[uG] = uM
		self.map_GM[vG] = vM
		self.map_MG[uM] = uG
		self.map_MG[vM] = vG
		

		self.edge_count[uG] += 1
		self.edge_count[vG] += 1
	
	def _revert_mappings(self, edge_idx):

		uG, vG, _ = self.G_edges[edge_idx]
		uM = self.map_GM.get(uG, -1)
		vM = self.map_GM.get(vG, -1)
		

		if uG in self.edge_count:
			self.edge_count[uG] -= 1
			if self.edge_count[uG] == 0:
				del self.edge_count[uG]
				self.map_GM[uG] = -1
				if uM != -1:
					self.map_MG[uM] = -1
		
		if vG in self.edge_count:
			self.edge_count[vG] -= 1
			if self.edge_count[vG] == 0:
				del self.edge_count[vG]
				self.map_GM[vG] = -1
				if vM != -1:
					self.map_MG[vM] = -1
	
	def count(self):

		G_edges = self.G_edges
		n_edges = len(self.M_edges)
		e_stack = self.e_stack
		edge_counter = self.edge_counter
		
		while True:
			next_index = self._find_next_match()
			self.eG = next_index
			
			if next_index < self.n:  
				if self.eM == n_edges - 1:  
					for idx in e_stack:
						u, v, _ = G_edges[idx]
						edge_counter[(u, v)] += 1
					u, v, _ = G_edges[next_index]
					edge_counter[(u, v)] += 1
	
					self.eG = next_index + 1
				else:  
					if not e_stack:
						_, _, t = G_edges[next_index]
						self.t0 = t + self.delta
					

					self._update_mappings(next_index)
					

					e_stack.append(next_index)
					self.eM += 1
					self.eG = next_index + 1
			else:  
				while self.eG >= self.n or (self.t0 is not None and self.eG < self.n and G_edges[self.eG][2] > self.t0):
					if not e_stack:
						return edge_counter
					pop_index = e_stack.pop()
					self.eG = pop_index + 1

					self._revert_mappings(pop_index)
					self.eM -= 1
					if not e_stack:
						self.t0 = None







