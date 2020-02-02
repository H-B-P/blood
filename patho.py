from fenics import *
# Create mesh and define function space
mesh = BoxMesh(Point(0, 0, 0),  Point(10, 10, 10), 10, 10, 10)
V = FunctionSpace(mesh, 'P', 1)

def boundary(x, on_boundary):
    return on_boundary

class WaveSurface(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary and near(x[2], 10):
            return True
        else:
            return False

class InputSurface(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary and near(x[0], 0) and between(x[1], (4, 6)) and between(x[2], (4, 6)):
            return True
        else:
            return False

class OutputSurface(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary and near(x[0], 10) and between(x[1], (4, 6)) and between(x[2], (4, 6)):
            return True
        else:
            return False

waveSurface=WaveSurface()
inputSurface=InputSurface()
outputSurface=OutputSurface()

boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
waveSurface.mark(boundaries, 1)
inputSurface.mark(boundaries, 2)
outputSurface.mark(boundaries, 3)

################################

meshFile = File("mesh.pvd")
meshFile<<mesh

boundariesFile = File("boundaries.pvd")
boundariesFile<<boundaries

################################

dt=0.1
g=9.81
rho=1000

#################


phi = TrialFunction(V)

prevPhi = Constant(0)
prevPhiVel = Constant(0)
prevPhiAcc = Constant(0)

phiVel = (phi-prevPhi)*Constant(2/dt)-prevPhiVel
phiAcc = (phiVel-prevPhiVel)*Constant(2/dt)-prevPhiAcc

v = TestFunction(V)

inFlow=Constant(1.0)
outFlow=Constant(0.0)
F = dot(grad(phi),grad(v))*dx

a,L = lhs(F),rhs(F)

phi = Function(V)
solve(a==L, phi)

phiFile = File("phi.pvd")
phiFile<<phi

###################

phiVel = (phi-prevPhi)*Constant(2/dt)-prevPhiVel
phiAcc = (phiVel-prevPhiVel)*Constant(2/dt)-prevPhiAcc

