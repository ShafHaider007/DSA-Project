"""Shortest Path Calculation using Dijkstra's Algorithm"""

This code calculates the shortest path between two points on a road network using Dijkstra's algorithm. The road network data is loaded from a GeoJSON file containing information about the road segments and their connectivity.

"""Dependencies"""

The code requires the following dependencies:

NetworkX: A Python library for the creation, manipulation, and study of the structure, dynamics, and functions of complex networks.
GeoPandas: A Python library for working with geospatial data, built on top of the pandas library.
Shapely: A Python library for manipulation and analysis of geometric objects.
heapq: A Python library providing an implementation of the heap queue algorithm.
json: A Python library for working with JSON data.
math: A Python library providing mathematical functions.

"""Functions"""

The code includes the following functions:

->haversine(lon1, lat1, lon2, lat2): Calculates the Haversine distance between two points given their longitudes and latitudes.

->dijkstra(graph, source, target): Applies Dijkstra's algorithm to find the shortest path in a graph.
Execution
To execute the code, follow these steps:
Provide the path to the GeoJSON file containing the road network data in the open() function where data is loaded.
Set the major_intersection_threshold variable to define the threshold for considering a node as a major intersection.
Define the source and destination nodes by specifying their longitude and latitude coordinates.

The code will calculate the shortest path between the source and destination nodes using Dijkstra's algorithm. The resulting shortest path nodes will be saved as a new GeoJSON file.

Note: If the source or destination nodes are not present in the road network data, the code will find the nearest node in the network and use it as the source or destination node.

""" Development"""
	 The Dijkstra algorithm get the shortest path and save the nodes in a geojson file containing the nodes that connecting with each other and making the shortest route.

Limitations

The code currently assigns edge weights based on a simple formula considering transportation mode and traffic conditions. 