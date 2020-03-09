from fenics import *
import math
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



gamma=0.5
dt=0.02

timespan=20
time=0

g=10
rho=1000



#################

c_one = (1-gamma)/(gamma*dt)
c_two = -1/((1-gamma)*gamma*dt)
c_three = (2-gamma)/((1-gamma)*dt)

#################



prevPhiFunc = Function(V)
prevPhiFuncVel = Function(V)
prevPhiFuncAcc = Function(V)

interPhiFunc = Function(V)
interPhiFuncVel = Function(V)
interPhiFuncAcc = Function(V)

phiFunc = Function(V)
phiFuncVel = Function(V)
phiFuncAcc = Function(V)



v = TestFunction(V)
phi = TrialFunction(V)

phiVel_first = (phi-prevPhiFunc)*2/(gamma*dt)-prevPhiFuncVel
phiAcc_first = (phiVel_first-prevPhiFuncVel)*2/(gamma*dt)-prevPhiFuncAcc

phiVel_second = c_one*prevPhiFunc + c_two*interPhiFunc + c_three*phi
phiAcc_second = c_one*prevPhiFuncVel + c_two*interPhiFuncVel + c_three*phiVel_second


pressure = Function(V)


inFlow=Constant(0)
outFlow=Constant(0)

F_first = dot(grad(phi),grad(v))*dx(0) + (1/float(g))*phiAcc_first*v*ds(1) - inFlow*v*ds(2) + outFlow*v*ds(3)
a_first,L_first = lhs(F_first),rhs(F_first)

F_second = dot(grad(phi),grad(v))*dx(0) + (1/float(g))*phiAcc_second*v*ds(1) - inFlow*v*ds(2) + outFlow*v*ds(3)
a_second,L_second = lhs(F_second),rhs(F_second)



phiFile = File("phi.pvd")

pressureFile = File("pressure.pvd")

while (time<timespan):
 
 inFlow.assign(Constant(math.sin(time)))
 
 #if (time<=1):
 # inFlow.assign(Constant(0))
 #else:
 # inFlow.assign(Constant(1))
 
 #step one
 
 solve(a_first==L_first, interPhiFunc)
 
 interPhiFuncVel.assign(((interPhiFunc-prevPhiFunc)*2/(dt*gamma))-prevPhiFuncVel)
 interPhiFuncAcc.assign(((interPhiFuncVel-prevPhiFuncVel)*2/(dt*gamma))-prevPhiFuncAcc)
 
 #step two
 
 solve(a_second==L_second, phiFunc)
 
 phiFuncVel.assign(c_one*prevPhiFunc + c_two*interPhiFunc + c_three*phiFunc)
 phiFuncAcc.assign(c_one*prevPhiFuncVel + c_two*interPhiFuncVel + c_three*phiFuncVel)
 
 #output output
 
 phiFile<<phiFunc
 
 pressure.assign(phiFuncVel*Constant(rho))
 pressureFile<<pressure
 
 print(round(time, 4), round(pressure([1,1,1]),2), round(pressure([1,1,2]),2))
 
 #prep for next step
 
 prevPhiFunc.assign(phiFunc)
 prevPhiFuncVel.assign(phiFuncVel)
 prevPhiFuncAcc.assign(phiFuncAcc)
 
 time+=dt
 
###################



