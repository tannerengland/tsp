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
		cities = self._scenario.getCities()

		doAgain = True
		while doAgain:
			random_index = random.randint(0, len(cities) - 1)
			route = []
			route.append(cities[random_index])
			currCity = cities[random_index]

			visitedBools = [False] * len(cities)
			visitedBools[random_index] = True

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
				currCity = cities[nextCityIndex]
				currIndex = nextCityIndex
				# set city to visited
				visitedBools[currIndex] = True

			if route[-1].costTo(route[0]) != float('inf'):
				doAgain = False

		return route



	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''

	def branchAndBound( self, time_allowance=60.0 ):
		pass



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

	def RMCA(self, edges, r, c, bound=0):
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
					if (row != r and col != c):
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
					if (row != c and col != r):
						bound = float('inf')
		ogMatrix = np.transpose(edgesTranspose)
		ogMatrixList = ogMatrix.tolist()
		# if fro < len(ogMatrix) and to < len(ogMatrix):
		# 	ogMatrix[fro][to] = float('inf')
		return bound, ogMatrixList

	values = [
		[float('inf'), 8, float('inf'), 4],
		[float('inf'), float('inf'), float('inf'), float('inf')],
		[2, 6, float('inf'), 4],
		[float('inf'), 3, float('inf'), float('inf')]
	]

	RMCA(0, values, 1, 2, 0)
