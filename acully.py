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

class BottomSurface(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary and near(x[2],0):
            return True
        else:
            return False


waveSurface=WaveSurface()
inputSurface=InputSurface()
outputSurface=OutputSurface()
bottomSurface=BottomSurface()

domains = CellFunction("size_t",mesh)
domains.set_all(0)

dx = Measure('dx', domain=mesh, subdomain_data=domains)

boundaries = FaceFunction("size_t", mesh)
boundaries.set_all(0)
waveSurface.mark(boundaries, 1)
inputSurface.mark(boundaries, 2)
outputSurface.mark(boundaries, 3)
bottomSurface.mark(boundaries,4)

ds = Measure("ds", domain=mesh, subdomain_data=boundaries)


################################

meshFile = File("mesh.pvd")
meshFile<<mesh

boundariesFile = File("boundaries.pvd")
boundariesFile<<boundaries

################################

dt=0.1
timespan=1

g=9.81
rho=1000

time=0

#################

prevPhi = Constant(0)
prevPhiVel = Constant(0)
prevPhiAcc = Constant(0)

phiFile = File("phi.pvd")

while (time<timespan):
 
 time+=dt
 
 phi = TrialFunction(V)
 
 phiVel = (phi-prevPhi)*Constant(2/dt)-prevPhiVel
 phiAcc = (phiVel-prevPhiVel)*Constant(2/dt)-prevPhiAcc
 
 v = TestFunction(V)
 
 if (time>0.5):
  inFlow=Constant(1.0)
 else:
  inFlow=Constant(0.0)
 
 outFlow=Constant(0.0)
 F = dot(grad(phi),grad(v))*dx(0) - g*phiAcc*v*ds(1) - inFlow*v*ds(2) + outFlow*v*ds(3)
 
 a,L = lhs(F),rhs(F)
 
 bc1 = DirichletBC(V, 0, boundaries, 4)
 bc2 = DirichletBC(V, 0, boundaries, 0)
 
 bcs = [bc1]#,bc2]
 
 phi = Function(V)
 solve(a==L, phi, bcs)
 
 phiFile<<phi
 
 phiVel = (phi-prevPhi)*Constant(2/dt)-prevPhiVel
 phiAcc = (phiVel-prevPhiVel)*Constant(2/dt)-prevPhiAcc
 
 prevPhi=phi
 prevPhiVel=phiVel
 prevPhiAcc=phiAcc

 
###################



