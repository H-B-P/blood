from fenics import *
from dolfin import *
from mshr import *
import math
# Create mesh and define function space

domain = Box(Point(0,0,0),Point(0.3,0.2,0.25)) + Cylinder(Point(-0.001,0.1,0.125),Point(0,0.1,0.125), 0.005, 0.005)
mesh = generate_mesh(domain,40)
V = FunctionSpace(mesh, 'P', 1)

def boundary(x, on_boundary):
    return on_boundary

class WaveSurface(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary and near(x[2], 0.25):
            return True
        else:
            return False

class InputSurface(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary and near(x[0], -0.001) and between(x[1], (0, 1)) and between(x[2], (0, 1)):
            return True
        else:
            return False

waveSurface=WaveSurface()
inputSurface=InputSurface()

domains = MeshFunction("size_t",mesh,3)
domains.set_all(0)

dx = Measure('dx', domain=mesh, subdomain_data=domains)

boundaries = MeshFunction("size_t", mesh,2)
boundaries.set_all(0)
waveSurface.mark(boundaries, 1)
inputSurface.mark(boundaries, 2)

ds = Measure("ds", domain=mesh, subdomain_data=boundaries)


################################

meshFile = File("mesh.pvd")
meshFile<<mesh

boundariesFile = File("boundaries.pvd")
boundariesFile<<boundaries

################################


gamma=0.5
dt=0.0000001

timespan=0.0001
time=0

g=9.81
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

F_first = dot(grad(phi),grad(v))*dx(0) + (1/float(g))*phiAcc_first*v*ds(1) - inFlow*v*ds(2)
a_first,L_first = lhs(F_first),rhs(F_first)

F_second = dot(grad(phi),grad(v))*dx(0) + (1/float(g))*phiAcc_second*v*ds(1) - inFlow*v*ds(2)
a_second,L_second = lhs(F_second),rhs(F_second)



phiFile = File("phi.pvd")

pressureFile = File("pressure.pvd")

record=[]

fives=[]
fifteens=[]
twentyfives=[]
fifties = [] 

while (time<(timespan+dt/2)):
 
 flowAmt = 4.5*math.cos(time*2*math.pi*95000)#0.05*math.cos(time*2*math.pi*1000000)+0.05*math.cos(time*2*math.pi*700000)+0.05*math.cos(time*2*math.pi*300000)
 print("flowAmt: ",flowAmt)
 inFlow.assign(Constant(flowAmt))
 
 #step one
 
 solve(a_first==L_first, interPhiFunc)
 
 interPhiFuncVel.assign(((interPhiFunc-prevPhiFunc)*2/(dt*gamma))-prevPhiFuncVel)
 interPhiFuncAcc.assign(((interPhiFuncVel-prevPhiFuncVel)*2/(dt*gamma))-prevPhiFuncAcc)
 
 #step two
 
 solve(a_second==L_second, phiFunc)
 
 phiFuncVel.assign(c_one*prevPhiFunc + c_two*interPhiFunc + c_three*phiFunc)
 phiFuncAcc.assign(c_one*prevPhiFuncVel + c_two*interPhiFuncVel + c_three*phiFuncVel)
 
 #output output
 
 step=round(time/dt)
 
 if (time>0):
  phiFile<<phiFunc
  
  pressure.assign(phiFuncVel*Constant(rho))
  pressureFile<<pressure
  
  addl=[pressure([0.003,0.1,0.125]), pressure([0.005,0.1,0.125]), pressure([0.007,0.1,0.125]), pressure([0.009,0.1,0.125]), pressure([0.011,0.1,0.125]), pressure([0.013,0.1,0.125]), pressure([0.015,0.1,0.125]), pressure([0.017,0.1,0.125]), pressure([0.019,0.1,0.125]), pressure([0.021,0.1,0.125]), pressure([0.023,0.1,0.125]), pressure([0.025,0.1,0.125]), pressure([0.027,0.1,0.125]), pressure([0.029,0.1,0.125])]
  
  record.append(addl)
  
  fives.append(round(pressure([0.005,0.1,0.125])))
  fifteens.append(round(pressure([0.015,0.1,0.125])))
  twentyfives.append(round(pressure([0.025,0.1,0.125])))
  fifties.append(round(pressure([0.05,0.1,0.125])))
  
 print(step, round(pressure([0.005,0.1,0.125]))/1000, "kPa", round(pressure([0.015,0.1,0.125]))/1000, "kPa", round(pressure([0.025,0.1,0.125]))/1000, "kPa", round(pressure([0.05,0.1,0.125]))/1000, "kPa")
  
 #prep for next step
 
 prevPhiFunc.assign(phiFunc)
 prevPhiFuncVel.assign(phiFuncVel)
 prevPhiFuncAcc.assign(phiFuncAcc)
 
 time+=dt


print(fives)
print(fifteens)
print(twentyfives)
print(fifties)
###################



