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

dt=0.01
timespan=10

g=10
rho=1000

time=0

#################

prevPhiFunc = Function(V)
prevPhiFuncVel = Function(V)
prevPhiFuncAcc = Function(V)

#prevPhiFunc.assign(Expression("0", degree=0))
#prevPhiFuncVel.assign(Expression("(10-x[2])*9.81", degree=1))
#prevPhiFuncAcc.assign(Expression("0", degree=0))

v = TestFunction(V)

phi = TrialFunction(V)
phiVel = ((phi-prevPhiFunc)*2/dt)-prevPhiFuncVel
phiAcc = (phiVel-prevPhiFuncVel)*Constant(2/dt)-prevPhiFuncAcc

phiFunc = Function(V)
phiFuncVel = Function(V)
phiFuncAcc = Function(V)

pressure = Function(V)

inFlow=Constant(0)
outFlow=Constant(0)

F = dot(grad(phi),grad(v))*dx(0) + (1/float(g))*phiAcc*v*ds(1) - inFlow*v*ds(2) + outFlow*v*ds(3)
a,L = lhs(F),rhs(F)

#bc1 = DirichletBC(V, 0, boundaries, 4)
#bc2 = DirichletBC(V, 0, boundaries, 0)
#bcs = [bc1,bc2]

solve(a==L, phiFunc)

phiFile = File("phi.pvd")

pressureFile = File("pressure.pvd")

while (time<timespan):
 
 #inFlow.assign(Constant(10*min(time, 0.1))) # This works after the ramp ends. (no it doesnt, numerical errors just cancel)
 #inFlow.assign(Constant(10*min(max(time-0.05, 0),0.1))) #This also only works under the surface once we hit 0.0
 
 inFlow.assign(Constant(1))
 
 solve(a==L, phiFunc)
 
 phiFile<<phiFunc
 
 phiFuncVel.assign(((phiFunc-prevPhiFunc)*2/dt)-prevPhiFuncVel)
 phiFuncAcc.assign(((phiFuncVel-prevPhiFuncVel)*2/dt)-prevPhiFuncAcc)
 
 pressure.assign(phiFuncVel*Constant(rho))
 
 print(round(time, 4), round(pressure([1,1,1]),2), round(pressure([1,1,2]),2))
 
 pressureFile<<pressure
 
 prevPhiFunc.assign(phiFunc)
 prevPhiFuncVel.assign(phiFuncVel)
 prevPhiFuncAcc.assign(phiFuncAcc)
 
 time+=dt

 
###################



