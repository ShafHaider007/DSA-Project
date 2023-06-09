import networkx as nx
import geopandas as gpd
import json
import math
from heapq import heappop, heappush
from shapely.geometry import Point
from shapely.ops import nearest_points

# Function to calculate Haversine distance between two points


def haversine(lon1, lat1, lon2, lat2):
    # Convert coordinates from degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * \
        math.cos(lat2) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    radius = 6371  # Radius of the Earth in kilometers
    distance = radius * c

    return distance

# Function to calculate the weight of an edge based on transportation mode and traffic information


def calculate_edge_weight(distance, transportation_mode, traffic_data):
    if transportation_mode == 'car':
        if traffic_data == 'congested':
            return distance * 1.5
        else:
            return distance
    elif transportation_mode == 'bike':
        return distance * 0.8
    elif transportation_mode == 'pedestrian':
        return distance * 0.7
    else:
        return distance

# Function to apply Dijkstra's algorithm for finding the shortest path


def dijkstra(graph, source, target):
    distances = {node: float('inf') for node in graph.nodes}
    distances[source] = 0
    previous = {node: None for node in graph.nodes}
    visited = set()
    heap = [(0, source)]

    while heap:
        current_distance, current_node = heappop(heap)

        if current_node == target:
            break

        if current_node in visited:
            continue

        visited.add(current_node)

        for neighbor in graph.neighbors(current_node):
            edge_weight = graph.edges[current_node, neighbor]['weight']
            distance = current_distance + edge_weight
            if distance < distances[neighbor]:
                distances[neighbor] = distance
                previous[neighbor] = current_node
                heappush(heap, (distance, neighbor))

    path = []
    current_node = target

    while current_node is not None:
        path.append(current_node)
        current_node = previous[current_node]

    path.reverse()

    return path


# Load the GeoJSON file
with open(r'C:\Users\shafs\OneDrive\Desktop\project\reduced_isb_roads.geojson', encoding='utf-8') as file:
    data = json.load(file)

# Create an empty graph
G = nx.Graph()

# Define a threshold for the number of incident edges to consider a node as a major intersection
major_intersection_threshold = 1

# Iterate over the features in the GeoJSON data
for feature in data['features']:
    geometry = feature['geometry']
    if geometry['type'] == 'LineString':
        coordinates = geometry['coordinates']
    elif geometry['type'] == 'MultiLineString':
        coordinates = [coord for sublist in geometry['coordinates']
                       for coord in sublist]
    else:
        continue

    # Iterate over the coordinate pairs and add edges to the graph
    for i in range(len(coordinates) - 1):
        u = tuple(coordinates[i])
        v = tuple(coordinates[i + 1])

        # Calculate the distance between u and v using the Haversine formula
        lon1, lat1 = u
        lon2, lat2 = v
        distance = haversine(lon1, lat1, lon2, lat2)

        # Assign weights to edges based on transportation mode and traffic conditions
        edge_weight = calculate_edge_weight(distance, 'car', 'congested')

        # Add the nodes and edges only if they satisfy the major intersection criteria
        if G.has_node(u):
            G.nodes[u]['incident_edges'] += 1
        else:
            G.add_node(u, incident_edges=1)
        if G.has_node(v):
            G.nodes[v]['incident_edges'] += 1
        else:
            G.add_node(v, incident_edges=1)
        if G.nodes[u]['incident_edges'] >= major_intersection_threshold:
            G.add_edge(u, v, weight=edge_weight)
        if G.nodes[v]['incident_edges'] >= major_intersection_threshold:
            G.add_edge(v, u, weight=edge_weight)

# Create a GeoPandas DataFrame from the GeoJSON data
gdf = gpd.GeoDataFrame.from_features(data)

# Calculate the bounding box from the GeoPandas DataFrame
xmin, ymin, xmax, ymax = gdf.total_bounds

# Define the source and destination nodes
source_node = (72.98863101900939, 33.64444206825158)
destination_node = (72.97951625020546, 33.645927458607865)

# Check if the source node is in the graph
if source_node not in G.nodes:
    # Find the nearest node in the graph to the given source node
    source_point = Point(source_node)
    nearest_node = min(
        G.nodes, key=lambda node: source_point.distance(Point(node)))
    source_node = nearest_node

# Check if the destination node is in the graph
if destination_node not in G.nodes:
    # Find the nearest node in the graph to the given destination node
    destination_point = Point(destination_node)
    nearest_node = min(
        G.nodes, key=lambda node: destination_point.distance(Point(node)))
    destination_node = nearest_node

# Apply Dijkstra's algorithm to find the shortest path
shortest_path = dijkstra(G, source_node, destination_node)

# Create a new GeoJSON dictionary for the shortest path nodes
shortest_path_geojson = {
    'type': 'FeatureCollection',
    'features': []
}

# Iterate over the nodes along the shortest path and create GeoJSON features
for node in shortest_path:
    feature = {
        'type': 'Feature',
        'geometry': {
            'type': 'Point',
            'coordinates': [node[0], node[1]]
        },
        'properties': {}
    }
    shortest_path_geojson['features'].append(feature)

# Save the shortest path nodes as a GeoJSON file
with open('shortest_path_nodes11.geojson', 'w') as file:
    json.dump(shortest_path_geojson, file)

# Print the shortest path nodes GeoJSON
# print(shortest_path_geojson)
