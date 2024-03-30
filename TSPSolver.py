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
import queue


# Partial Path object to store cost matrix, cost, and cities that have been visited
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

    # assures that the priority queue prioritizes based off of the minimum of both cost and amount of cities visited
    def __lt__(self, other):
        return self._cost/len(self._visited_cities) < other._cost/len(other._visited_cities)


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None

    def setupWithScenario(self, scenario):
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

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time() - start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
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

    # A greedy approach to the traveling salesperson problem by finding the minimum outgoing edge to each unvisited city
    # until all cities are visited
    # Time: O(n^2) we must go through all the cities and get their associated costs to every city to decided which city
    # to visit next
    # Space: O(n) we must store the cities which is of length n (for the original cities array, the visited boolean
    # array, the route, and the costs of all the cities)
    def greedy(self, time_allowance=60.0):
        cities = self._scenario.getCities().copy()
        results = {}
        bssf = None

        start_time = time.time()
        doAgain = True
        count = 0
        while doAgain:
            count += 1
            # starts at a random city
            random_index = random.randint(0, len(cities) - 1)
            route = []
            route.append(cities[random_index])
            currCity = cities[random_index]

            # marks all cities as unvisited except the one we are starting from
            visitedBools = [False] * len(cities)
            visitedBools[random_index] = True
            cost = 0

            # while all cities are not visited
            while not all(visitedBools):
                costs = []
                i = 0
                for city in cities:
                    # if city is not visited, store cost
                    if not visitedBools[i]:
                        costs.append(currCity.costTo(city))
                    else:
                        costs.append(float('inf'))
                    i += 1
                # finds city that has the minimum cost of all unvisited cities and adds it to our route
                nextCityIndex = np.argmin(costs)
                route.append(cities[nextCityIndex])
                cost += costs[nextCityIndex]
                # sets the minimum cost city as current city to iterate and repeat
                currCity = cities[nextCityIndex]
                currIndex = nextCityIndex
                costs[currIndex] = float('inf')
                # set city to visited
                visitedBools[currIndex] = True

            # assures we have a path back to the first city we visited
            if route[-1].costTo(route[0]) != float('inf') and cost != float('inf'):
                cost += route[-1].costTo(route[0])
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

    # Branch and Bound algorithm that finds the minimum path amongst cities by considering paths and pruning those
    # that are not better than our Best Solution we have calculated So Far
    # Time: O(n^2 * b) where b is the number of states that we have made and n is the number of cities. The time it
    # takes for the number of nodes created for the potential queue insertion multiplied by the complexity of creating a
    # node which has a cost of O(n^2) (from the reduced cost matrix).
    # Space: O(n^2 * b) where b is the number of states that we have made and n is the number of cities. Since we store
    # a reduced cost matrix which is an nxn matrix for every state we have created.
    def branchAndBound(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities().copy()
        edges = [[None] * len(cities) for _ in range(len(cities))]
        # creates initial matrix of costs from city to city
        for i in range(len(cities)):
            for j in range(len(cities)):
                edges[i][j] = cities[i].costTo(cities[j])

        start_time = time.time()

        numCities = len(cities)
        pruned = 0
        count = 0
        # reduces our initial matrix and finds bound
        bound, matrix = self.lowerbound(edges, numCities, numCities, 0, False)
        # marks first city as visited
        visitedCities = []
        visitedCities.append(0)
        visitedCitiesNoIndex = []
        visitedCitiesNoIndex.append(cities[0])
        bssf = TSPSolution(visitedCitiesNoIndex)
        # creates initial partial path object
        P0 = PartialPath(matrix, visitedCities, bound)
        total = 1
        # initializes priority queue and pushes inital path on
        S = []
        heapq.heappush(S, P0)
        max = 1
        # finds a rough bssf based off of my greedy algorithm
        bssf.cost = self.greedy()['cost']

        # while the queue is not empty and hasn't been computing for longer than time_allowance, perform branch & bound
        while S and (time_allowance > (time.time() - start_time)):
            P = heapq.heappop(S)

            # if current path cost is less than our best solution so far, find children paths
            if P.get_cost() < bssf.cost:
                T = self.expand(P, P.get_visited_cities()[len(P.get_visited_cities()) - 1])
                total += len(T)
                pruned += numCities - len(P.get_visited_cities()) - len(T)
                for Pi in T:
                    # for each child, see if it is a complete path and less than our current BSSF, if so, update BSSF
                    if self.test(Pi, cities) < bssf.cost:
                        bssf.cost = self.test(Pi, cities)
                        count += 1
                        visitedCitiesNoIndex = []
                        for city in Pi.get_visited_cities():
                            visitedCitiesNoIndex.append(cities[city])
                        route = visitedCitiesNoIndex
                        temp = TSPSolution(route)
                        bssf = temp
                    # if a child's cost is less than bssf but not a complete path, push to the queue for further
                    # investigation of path
                    elif (Pi.get_cost() < bssf.cost):
                        heapq.heappush(S, Pi)
                        if len(S) > max:
                            max = len(S)
                    # otherwise prune
                    else:
                        pruned += 1
        end_time = time.time()

        results['cost'] = bssf.cost
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = max
        results['total'] = total
        results['pruned'] = pruned

        return results


    # Finds all children matrices/costs (or visits all possible cities) and return it in a list
    # Time: O(n^3) must calculate and create an nxn matrix of all costs for each child of a parent via the lowerbound
    # function for n-1 times worst case
    # Space: O(n^3) creates reduced cost matrices for all children of a given parent, which is an nxn dimensioned table
    # via the lowerbound function for n-1 times worst case
    def expand(self, parent, currRow):
        childrenPaths = []

        for col in range(len((parent.get_matrix()[0]))):
            # makes sure we are not trying to visit an already visited city or our self
            if (parent.get_matrix()[currRow][col] != float('inf') and currRow != col) and (col not in parent.get_visited_cities()):
                currBound, currEdges = self.lowerbound(parent.get_matrix(), currRow, col, parent.get_cost(), True)
                visitedCities = copy.deepcopy(parent.get_visited_cities())
                # marks city as visited for child
                visitedCities.append(col)
                currChild = PartialPath(currEdges, visitedCities, currBound)
                # appends child to list to be iterated through
                childrenPaths.append(currChild)

        return childrenPaths

    # Updates the cost matrix and does our reduced cost matrix algorithm, returning the lowerbound and cost matrix
    # Time: O(n^2) must go through the entirety of the nxn matrix to calculate the lower bound and the reduced cost
    # matrix (the minimum of each row and column)
    # Space: O(n^2) we must make and modify a deep copy of the matrix we are given so that we don't edit the parent
    # matrix to make the child and update the lowerbound
    def lowerbound(self, edges, r, c, prevBound, updateMat):
        edges = copy.deepcopy(edges)

        # tells us whether we should update a matrix
        if updateMat == True:
            # this adds path costs
            prevBound += edges[r][c]
            # this infinities out
            for row in range(len(np.transpose(edges)[0])):
                edges[row][c] = float('inf')
            for col in range(len(edges[0])):
                edges[r][col] = float('inf')
            # infinities out inverse
            edges[c][r] = float('inf')

        # finds bound of/minimizes rows
        bound = 0
        for row in range(len(edges[0])):
            minimum = min(edges[row])
            if minimum != float('inf'):
                bound += minimum
                for col in range(len(edges[0])):
                    edges[row][col] -= minimum

        # finds bound of/minimizes columns
        edgesTranspose = np.transpose(edges)
        for row in range(len(edgesTranspose[0])):
            minimum = min(edgesTranspose[row])
            if minimum != float('inf'):
                bound += minimum
                for col in range(len(edgesTranspose[0])):
                    edgesTranspose[row][col] -= minimum
        ogMatrix = np.transpose(edgesTranspose)
        ogMatrixList = ogMatrix.tolist()

        # adds previous bound to our current bound
        bound += prevBound

        return bound, ogMatrixList

    # Returns a partial path's cost if all cities are visited for that path
    # Time: O(n) we must iterate through the entirety of the amount of cities at worst case to check if all are visited
    # Space: O(1) we are returning the cost retrieved from a partial path object
    def test(self, path, cities):
        # makes sure all cities are visited to return the partial path cost
        if all(i in path.get_visited_cities() for i in range(0, len(cities))):
            return path.get_cost()
        else:
            return float('inf')

    ''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

    def fancy(self, time_allowance=60.0):
        pass