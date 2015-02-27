import matplotlib.pyplot as plt
def plot(shapelyGeometries):
    'plot shapelyGeometries'
    figure = plt.figure(num=None, figsize=(4, 4), dpi=180)
    axes = plt.axes()
    axes.set_aspect('equal', 'datalim')
    axes.xaxis.set_visible(False)
    axes.yaxis.set_visible(False)
    draw(shapelyGeometries)
            
def draw(gs):
    'Draw shapelyGeometries'
    # Handle single and multiplte geometries
    try:
        gs = iter(gs)
    except TypeError:
        gs = [gs]
    # For each shapelyGeometry,
    for g in gs:
        gType = g.geom_type
        if gType.startswith('Multi') or gType == 'GeometryCollection':
            draw(g.geoms)
        else:
            draw_(g)
            
def draw_(g):
    'Draw a shapelyGeometry; thanks to Sean Gilles'
    gType = g.geom_type
    if gType == 'Point':
        plt.pltot(g.x, g.y, 'k,')
    elif gType == 'LineString':
        x, y = g.xy
        plt.plot(x, y, 'b-')
    elif gType == 'Polygon':
        x, y = g.exterior.xy
        plt.fill(x, y, color='#cccccc', aa=True) 
        plt.plot(x, y, color='#666666', aa=True, lw=1.0)
        for hole in g.interiors:
            x, y = hole.xy
            plt.fill(x, y, color='#ffffff', aa=True) 
            plt.plot(x, y, color='#999999', aa=True, lw=1.0)