from fenics import *
# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)
# Define boundary condition
u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)
#
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx




# Compute solution
u = Function(V)
solve(a == L, u, bc)
# Plot solution and mesh
#plot(u)
#plot(mesh)
# Save solution to file in VTK format
vtkfile = File("poisson/solution.pvd")
vtkfile << u
