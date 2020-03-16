from dolfin import *
from fenics import *
from mshr import *
domain = Rectangle(Point(0., 0.), Point(5., 5.)) - Rectangle(Point(2., 1.25), Point(3., 1.75)) - Circle(Point(1, 4), .25) - Circle(Point(4, 4), .25)
domain.set_subdomain(1, Rectangle(Point(1., 1.), Point(4., 3.)))
domain.set_subdomain(2, Rectangle(Point(2., 2.), Point(3., 4.)))

mesh = generate_mesh(domain, 50)

meshFile = File("mesh3.pvd")
meshFile<<mesh
