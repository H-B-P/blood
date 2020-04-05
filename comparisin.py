import plotly.graph_objects as go
import numpy as np
import math



fig = go.Figure(
    data=[go.Line(x=list(np.arange(0.2,100, 0.2)),y=2.1*np.sin(math.pi*2*0.095*np.arange(0.2,100,0.2)))]
)
fig.show()
