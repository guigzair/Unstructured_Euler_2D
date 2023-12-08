import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from matplotlib import tri as mtri
from mpl_toolkits.axes_grid1 import make_axes_locatable
import meshpy.triangle as triangle

# Parameters
gamma = 1.4  # Specific heat ratio
dt = 0.0001   # Time step
t_end = 0.01  # End time
Lx, Ly = 1.0, 1.0  # Domain size

N = 250
##############################################################################
##############    Mesh #######################################################
##############################################################################
# Generate a random set of points representing the unstructured mesh
# np.random.seed(42)
num_points = 100
points = np.column_stack((np.random.rand(num_points) * Lx, np.random.rand(num_points) * Ly))
points = np.append(points, [[0,0],[Lx,0],[0,Ly],[Lx,Ly]], axis = 0)
# Create a Delaunay triangulation
triangulation = Delaunay(points)

def round_trip_connect(start, end):
    result = []
    for i in range(start, end):
        result.append((i, i+1))
    result.append((end, start))
    return result

boundaries = [[0,0],[Lx,0],[Lx,Ly],[0,Ly]]
info = triangle.MeshInfo()
info.set_points(boundaries)
info.set_facets(round_trip_connect(0, len(boundaries)-1))
mesh = triangle.build(info, max_volume=3e-2, min_angle=25)
mesh_points = np.array(mesh.points)
mesh_tris = np.array(mesh.elements)
mesh_facets = np.array(mesh.facets)
mesh_tris_neighbors = np.array(mesh.neighbors)

def getCenterTriangles(mesh_tris, mesh_points):
	barycenter = np.array([np.mean(mesh_points[mesh_tris[i]], 0) for i in range(len(mesh_tris))])
	return barycenter

def getNormals(mesh_tris, mesh_points):
	normals = np.zeros((len(mesh_tris),3,2))
	for tri in range(len(mesh_tris)):
		points_tri = mesh_tris[tri]
		for j in range(3):
			nx, ny = mesh_points[points_tri[(j + 1) % 3]] - mesh_points[points_tri[j]]
			normals[tri,j,0], normals[tri,j,1] = ny, -nx  # Rotate 90 degrees (assuming counterclockwise node order)
	return normals

mesh_baryTri = getCenterTriangles(mesh_tris, mesh_points)

mesh_normals = getNormals(mesh_tris, mesh_points)

##############################################################################
##############    Solving   ##################################################
##############################################################################

# Initial conditions
primitives = np.zeros((4,len(mesh_tris)))
for i in range(len(mesh_tris)):
	if np.linalg.norm(mesh_baryTri[i] - np.array([0.5,0.5]), 2) < 0.1:
		prim = np.array([1.,0.,0.,1.])
	else : 
		prim = np.array([0.125,0.,0.,0.1])
	primitives[:,i] = prim

# primitives = np.array([rho,u,v,p])

def getFlux(primitives_L, primitives_R, nx, ny, gamma = 1.4):
	# Using a simple upwind scheme for advection terms
	rho_L = primitives_L[0]
	u_L = primitives_L[1]
	v_L = primitives_L[2]
	P_L = primitives_L[3]

	rho_R = primitives_R[0]
	u_R = primitives_R[1]
	v_R = primitives_R[2]
	P_R = primitives_R[3]

    # Maximum wavelenghts
	C_L = np.sqrt(np.abs(gamma*P_L/np.maximum(rho_L,0.001))) + np.abs(u_L)
	C_R = np.sqrt(np.abs(gamma*P_L/np.maximum(rho_L,0.001))) + np.abs(u_L)
	C = np.maximum( C_L, C_R )
	
	# Energy
	en_L = P_L/(gamma-1)+0.5*rho_L * (u_L**2 + v_L**2)
	en_R = P_R/(gamma-1)+0.5*rho_R * (u_R**2 + v_R**2)

	# Flux
	flux_rho_L = rho_L * (u_L * nx + v_L * ny)
	flux_ru_L = rho_L * ( u_L**2 * nx + u_L*v_L * ny) + P_L * nx
	flux_rv_L = rho_L * (u_L*v_L * nx + v_L**2 * ny) + P_L * ny
	flux_E_L = (en_L + P_L) * (u_L * nx + v_L * ny)

	flux_rho_R = rho_R * (u_R * nx + v_R * ny)
	flux_ru_R = rho_R * ( u_R**2 * nx + u_R*v_R * ny) + P_R * nx
	flux_rv_R = rho_R * (u_R*v_R * nx + v_R**2 * ny) + P_R * ny
	flux_E_R = (en_R + P_R) * (u_R * nx + v_R * ny)

	# Total flux
	flux_rho = (flux_rho_L + flux_rho_R)/2 - C * 0.5 * (rho_R - rho_L)
	flux_ru =(flux_ru_L + flux_ru_R)/2 - C * 0.5 * (rho_R * u_R - rho_L * u_L)
	flux_rv = (flux_rv_L + flux_rv_R)/2 - C * 0.5 * (rho_R * v_R - rho_L * v_L)
	flux_E = (flux_E_L + flux_E_R)/2 - C * 0.5 * (en_R - en_L)

	return np.array([flux_rho, flux_ru, flux_rv, flux_E])
	

def getConserved(Primitives, gamma = 1.4):
	rho = Primitives[0,...]
	u = Primitives[1,...]
	v = Primitives[2,...]
	P = Primitives[3,...]
	Mass   = rho
	Mom_x   = rho * u 
	Mom_y   = rho * v 
	Energy = (P/(gamma-1) + 0.5*rho*(u**2+v**2))
	W = np.stack([Mass, Mom_x, Mom_y, Energy], axis = 0)
	return W


def getPrimitive(W, gamma = 1.4):
	Mass = W[0,...]
	Mom_x = W[1,...]
	Mom_y = W[2,...]
	Energy = W[3,...]
	rho = Mass
	rho = np.maximum(rho, 0.0001)
	u  = Mom_x / rho 
	v  = Mom_y / rho 
	P   = (Energy - 0.5*rho * (u**2 + v**2)) * (gamma-1)
	Primitives = np.stack([rho, u, v, P],axis = 0)
	return Primitives

# Prim = np.array([rho, u, v, p])
# W = getConserved(np.array([rho, u, v, p]))
# primitives = getPrimitive(W)

# Time-stepping loop
t = 0
for _ in range(N):
    # Compute fluxes (you would need a proper numerical flux function)
	rho = primitives[0,...]
	u = primitives[1,...]
	v = primitives[2,...]
	P = primitives[3,...]

	W = getConserved(primitives)
	W_updated = np.zeros(W.shape)
    # Update solution using a simple upwind scheme
	for tri in range(len(mesh_tris)):
		i, j, k = mesh_tris[tri]
		neighbors = mesh_tris_neighbors[tri]

		area = 0.5 * ((mesh_points[j, 1] - mesh_points[i, 1]) * (mesh_points[k, 0] - mesh_points[i, 0]) - 
                      (mesh_points[j, 0] - mesh_points[i, 0]) * (mesh_points[k, 1] - mesh_points[i, 1]))

		diff = 0
		# For loop around edges of the triangle
		for n in range(3):
			# Get edge
			pt_0 = mesh_points[mesh_tris[tri, n]]
			pt_1 = mesh_points[mesh_tris[tri, (n+1)%3]]

			# Find neighbor who has the edge in common
			neigh = -1
			for ne in neighbors:
				if len(np.intersect1d(np.array([mesh_tris[tri, n], mesh_tris[tri, (n+1)%3]]),mesh_tris[ne])) == 2:
					neigh = ne
			
			# Boundary condition
			if neigh == -1:
				neigh_primitives = primitives[:,tri]
			else:
				neigh_primitives = primitives[:,neigh]
			
			# Find normal of edge
			nx, ny = pt_1 - pt_0
			nx, ny = ny, -nx

			# Compute flux
			F = getFlux(primitives[:,tri], neigh_primitives, nx, ny, gamma = 1.4)
			diff += F

		W_updated[:,tri] = W[:,tri] - dt / np.maximum(np.abs(area),1e-10) * diff
	primitives = getPrimitive(W_updated)

	t += dt

##############################################################################
##############    Plotting   #################################################
##############################################################################

# Plot the results
triang = mtri.Triangulation(mesh_points[:, 0], mesh_points[:, 1], mesh_tris)

# fig, ax = plt.subplots(1, 1, figsize=(12, 8))
# ax.cla()
# ax.set_aspect('equal')
# tpc = ax.tricontourf(triang, rho, cmap='viridis')
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.05)
# fig.colorbar(tpc,cax = cax)

# ax.triplot(triang, 'ko-', ms=0.5, lw=0.3)
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_title('Density')



fig, ax = plt.subplots()
ax.set_aspect('equal')
tpc = ax.tripcolor(triang, facecolors=primitives[0], edgecolors='k')
fig.colorbar(tpc)
plt.show()


