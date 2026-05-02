from pylab import *

a = 1.1
steps = 30

def init():
    global x, result
    x = 0.1
    result = [x]
    
def observe():
    global x, result
    result.append(x)
    
def f(x): # iterative map
    return x + 2.5*x * (1-x)

def update():
    global x, result
    x = f(x)
    
init()
for t in range(steps):
    update()
    observe()
    
xmin, xmax = 0, 1.5
plot([xmin, xmax], [xmin, xmax], 'k')

rng = arange(xmin, xmax, (xmax-xmin)/100.)
plot(rng, f(rng), 'k')

horizontal = [result[0]]
vertical = [result[0]]

for x in result[1:]:
    horizontal.append(vertical[-1])
    vertical.append(x)
    horizontal.append(x)
    vertical.append(x)
plot(horizontal, vertical, 'b')

show()
