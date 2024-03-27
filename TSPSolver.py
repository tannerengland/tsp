#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))


import time
import numpy as np
from TSPClasses import *
import heapq
import itertools
import copy

class PartialPath:
    def __init__(self, matrix, visited_cities=None, cost=0):
        self._matrix = matrix
        self._visited_cities = visited_cities or []
        self._cost = cost

    def get_matrix(self):
        return self._matrix

    def get_visited_cities(self):
        return self._visited_cities

    def get_cost(self):
        return self._cost

class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution,
		time spent to find solution, number of permutations tried during search, the
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	def greedy( self,time_allowance=60.0 ):
		cities = self._scenario.getCities().copy()
		results = {}
		bssf = None

		start_time = time.time()
		doAgain = True
		count = 0
		while doAgain:
			count += 1
			random_index = random.randint(0, len(cities) - 1)
			route = []
			route.append(cities[random_index])
			currCity = cities[random_index]
			currIndex = random_index

			visitedBools = [False] * len(cities)
			visitedBools[random_index] = True
			cost = 0

			# while all cities are not visited
			while not all(visitedBools):
				# costs = np.array(len(cities))
				costs = []
				i = 0
				for city in cities:
					#if city is not visited
					if not visitedBools[i]:
						costs.append(currCity.costTo(city))
					else:
						costs.append(float('inf'))
					i += 1
				nextCityIndex = np.argmin(costs)
				route.append(cities[nextCityIndex])
				cost += costs[nextCityIndex]
				currCity = cities[nextCityIndex]
				currIndex = nextCityIndex
				costs[currIndex] = float('inf')
				# set city to visited
				visitedBools[currIndex] = True

			if route[-1].costTo(route[0]) != float('inf') and cost != float('inf'):
				doAgain = False

			bssf = TSPSolution(route)
			bssf.cost = cost

		end_time = time.time()
		results['cost'] = bssf.cost
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results



	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''

	def branchAndBound( self, time_allowance=60.0 ):
		cities = self._scenario.getCities().copy()
		# distances = self.cities[0].costTo(cities[1]._index)
		# dict = {}
		edges = [[None] * len(cities) for _ in range(len(cities))]

		for i in range(len(cities)):
			for j in range(len(cities)):
				edges[i][j] = cities[i].costTo(cities[j])
				# if j < i:
				# 	i,j = j,i
				# 	dict[(i, j)] = self.cities[i].costTo(cities[j]._index)

		# bssf, edges = self.partial_path_tree(edges)
# lowerbound is RCMA
		# expand is partial path
		# test sees if all cities are visited

		visitedCities = [False] * len(cities)
		visitedCities[0] = True
		bound, matrix = self.lowerbound(edges, len(edges[0]), len(edges[0]), 0, [], [])
		P0 = PartialPath(matrix, visitedCities, bound)
		S = [P0]
		bssf = self.greedy()['cost']
		while S:
			P = heapq.heappop(S)
			lowerbound, _ = self.lowerbound(P.get_matrix, len(P.get_matrix), len(P.get_matrix), P.get_cost, [], [])
			if lowerbound < bssf:
				T = self.expand(P)

		return bssf, edges

	# expand will be finding all children's matrices and costs and return that in a list (update city to visited here)
	# test will compare cost of a particular matrix that has visited Cities all equal to true (if the city is the last city in the array, then make sure there is a valid path to first city)

	def expand(self, parent, currRow):
		childrenPaths = []

		for col in range(len((parent.get_matrix[0]))):
			if (parent.get_matrix[currRow][col] != float('inf') and currRow != col):
				currBound, currEdges = self.partial_path(parent.get_matrix, currRow, col, parent.get_cost, prevRows, prevCols)
				visitedCities = parent.get_visitedCities()
				visitedCities[col] = True
				currChild = PartialPath(currEdges, visitedCities, currBound)
				childrenPaths.append(currChild)
			else:
				currBound = float('inf')
				currEdges = "place_holder"
				visitedCities = parent.get_visitedCities()
				# visitedCities[col] = True
				currChild = PartialPath(currEdges, visitedCities, currBound)
				childrenPaths.append(currChild)

		return childrenPaths

	def partial_path(self, mat, r, c, prevBound, prevRows, prevCols):
		edges = copy.deepcopy(mat)

		prevBound += edges[r][c]
		for row in range(len(np.transpose(edges)[0])):
			for col in range(len(edges[0])):
				if (row == r or col == c):
					edges[row][col] = float('inf')
		edges[c][r] = float('inf')

		bound, newEdges = self.lowerbound(edges, r, c, prevBound, prevRows, prevCols)
		return bound, newEdges

	def lowerbound(self, edges, r, c, prevBound, prevRows, prevCols):
		bound = 0
		for row in range(len(edges[0])):
			boundAdded = False
			minimum = min(edges[row])
			for col in range(len(edges[0])):
				# (minimum == float('inf') and (row != fro and col != to)) or
				if minimum != float('inf'):
					edges[row][col] -= minimum
					if boundAdded == False:
						bound += minimum
						boundAdded = True
				else:
					if (row != r and col != c) and (row not in prevRows and col not in prevCols):
							bound = float('inf')
		edgesTranspose = np.transpose(edges)
		for row in range(len(edgesTranspose[0])):
			boundAdded = False
			minimum = min(edgesTranspose[row])
			for col in range(len(edgesTranspose[0])):
				# (minimum == float('inf') and (row != fro and col != to)) or
				if minimum != float('inf'):
					edgesTranspose[row][col] -= minimum
					if boundAdded == False:
						bound += minimum
						boundAdded = True
				else:
					if (row != c and col != r) and (row not in prevCols and col not in prevRows):
						bound = float('inf')
		ogMatrix = np.transpose(edgesTranspose)
		ogMatrixList = ogMatrix.tolist()
		# if fro < len(ogMatrix) and to < len(ogMatrix):
		# 	ogMatrix[fro][to] = float('inf')
		bound += prevBound
		return bound, ogMatrixList

	def test(bool_list):
		return all(bool_list)

	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

	def fancy( self,time_allowance=60.0 ):
		pass

	# def partial_path_tree(self, ogMat):
	# 	edges = copy.deepcopy(ogMat)
	# 	done = False
	#
	# 	prevRows = []
	# 	prevCols = []
	# 	bssf, initMat = self.lowerbound(edges, len(edges), len(edges), 0, [], [])
	# 	edges = initMat
	# 	currRow = 0
	# 	counter = 0
	# 	bssfBefore = bssf
	# 	while not done:
	# 		bssfBefore = bssf
	# 		counter += 1
	#
	# 		potentialBounds = []
	# 		potentialEdges = []
	# 		for col in range(len((edges[0]))):
	# 			# edges[currRow][col] != float('inf') and
	# 			if (edges[currRow][col] != float('inf') and currRow != col):
	# 				currBound, currEdges = self.expand(edges, currRow, col, bssf, prevRows, prevCols)
	# 				potentialBounds.append(currBound)
	# 				potentialEdges.append(currEdges)
	# 			else:
	# 				potentialBounds.append(float('inf'))
	# 				potentialEdges.append("place_holder")
	# 		nextIndex = np.argmin(potentialBounds)
	# 		bssf = potentialBounds[nextIndex]
	# 		edges = potentialEdges[nextIndex]
	#
	# 		prevRows.append(currRow)
	# 		prevCols.append(nextIndex)
	#
	# 		currRow = nextIndex
	#
	# 		if bssfBefore < bssf:
	# 			done = True
	#
	# 		# TODO: make partial path class, make a cost matrix, cities that have already been visited (see if # of cities == # of cities) THEN put in PQ (already built)
	#
	# 	return bssfBefore, edges

# values = [
# 	[float('inf'), 1, float('inf'), 0, float('inf')],
# 	[float('inf'), float('inf'), 1, float('inf'), 0],
# 	[float('inf'), 0, float('inf'), 1, float('inf')],
# 	[float('inf'), 0, 0, float('inf'), 6],
# 	[0, float('inf'), float('inf'), 9 ,float('inf')]
# ]
# values = [
# 	[float('inf'), float('inf'), float('inf'), float('inf'), float('inf')],
# 	[float('inf'), float('inf'), 1, float('inf'), 0],
# 	[float('inf'), 0, float('inf'), float('inf'), float('inf')],
# 	[float('inf'), 0, 0, float('inf'), 6],
# 	[0, float('inf'), float('inf'), float('inf'), float('inf')]
# ]
# values = [
# 	[float('inf'), float('inf'), float('inf'), float('inf'), float('inf')],
# 	[float('inf'), float('inf'), float('inf'), float('inf'), 0],
# 	[float('inf'), 0, float('inf'), float('inf'), float('inf')],
# 	[float('inf'), float('inf'), float('inf'), float('inf'), float('inf')],
# 	[0, float('inf'), float('inf'), float('inf'), float('inf')]
# ]


values = [
	[float('inf'), 7, 3, 12],
	[3, float('inf'), 6, 14],
	[5, 8, float('inf'), 6],
	[9, 3, 5, float('inf')],
]
# values = [
# 	[float('inf'), 4, 0, 8],
# 	[0, float('inf'), 3, 10],
# 	[0, 3, float('inf'), 0],
# 	[6, 0, 2, float('inf')],
# ]
# values = [
# 	[float('inf'), float('inf'), float('inf'), float('inf')],
# 	[0, float('inf'), float('inf'), 10],
# 	[float('inf'), 3, float('inf'), 0],
# 	[6, 0, float('inf'), float('inf')],
# ]
# values = [
# 	[float('inf'), float('inf'), float('inf'), float('inf')],
# 	[float('inf'), float('inf'), float('inf'), float('inf')],
# 	[float('inf'), 3, float('inf'), 0],
# 	[float('inf'), 0, float('inf'), float('inf')],
# ]

# solver_object = TSPSolver(None)  # Replace value1 and value2 with actual parameters
# TSPSolver.partial_path_tree(solver_object, values)
# TSPSolver.partial_path(solver_object, values, 2, 3, 15, [0,1], [2,0])

		# RMCA(0, values, 1, 2, 0)
