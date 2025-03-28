import numpy as np

def grab_z(x, y, data):
    x = np.round(x, 2)
    y = np.round(y, 2)
    xs = np.where(data[:,0] == x)
    ys = np.where(data[:,1] == y)
    
    intersect = np.intersect1d(xs, ys)
    if len(intersect) == 0:
        print(x,y)
        return None
    else:
        return data[intersect[0], 2]

def gradient_LSQ(pointA, pointB, pointC, concentrationFunction):
    """Builds linear system to find gradient at point B
    """
    fa = grab_z(*pointA, concentrationFunction)
    fb = grab_z(*pointB, concentrationFunction)
    fc = grab_z(*pointC, concentrationFunction)
    xA, yA = pointA
    xB, yB = pointB
    xC, yC = pointC
    A = np.array([[xA - xB, yA - yB],
                  [xC - xB, yC - yB]])
    b = np.array([fa - fb, fc - fb])

    x = np.linalg.pinv(A) @ b

    return x

def hessian_LSQ(pointA, 
                pointB, 
                pointC, 
                pointD,
                pointE,
                pointF,
                pointG,
                pointH,
                pointI,
                gradientE, 
                concentrationFunction):
    """Builds linear system to find Hessian at point C, which we
    are assuming is the point of crossing the critical path. We can
    use the gradient we calculated at the point of crossing. 
    """
    xE = pointE[0]
    yE = pointE[1]
    fxe, fye = gradientE
    fe = grab_z(*pointE, concentrationFunction)

    A = np.zeros((8,3))
    b = np.zeros((8,1))
    for i, point in enumerate([pointA, pointB, pointC, pointD,
                               pointF, pointG, pointH, pointI]):
        x, y = point
        f = grab_z(x, y, concentrationFunction)
        A[i] = np.array([(x - xE)**2, (x-xE)*(y-yE), (y-yE)**2])
        b[i] = f - fe - fxe*(x-xE) - fye*(y-yE)

    x = np.linalg.pinv(A) @ b
    x = x.flatten()
    
    H = np.array([[x[0], x[1]], 
                  [x[1], x[2]]])

    return H

def unit_tangent_LSQ(gradB):
    tangent = np.array([-gradB[0], gradB[1]])
    norm = np.linalg.norm(tangent, ord=2)
    return tangent/norm


def calc_curvature(g, H):
    Fx = g[0]
    Fy = g[1]
    Fxx = H[0,0]
    Fxy = H[0,1]
    Fyy = H[1,1]
    num = -Fy**2*Fxx + 2*Fx*Fy*Fxy - Fx**2*Fyy
    denom = (Fx**2 + Fy**2)**(3/2)
    return num / denom

def calc_curvature_LSQ(pointA,
                       pointB,
                       pointC,
                       pointD,
                       pointE,
                       pointF,
                       pointG,
                       pointH,
                       pointI,
                       concentrationFunction):
    """Wrapper function to calculate curvature at point B
    """

    gradB = gradient_LSQ(pointA, pointB, pointC, concentrationFunction)
    H = hessian_LSQ(pointA, pointB, pointC,
                       pointD, pointE, pointF,
                       pointG, pointH, pointI,
                       gradB, concentrationFunction)
    curvature = calc_curvature(gradB, H)

    return gradB, curvature


def curvature_derivative(curvature_vals):
    """Calculate the derivative of the curvature
    """
    x = curvature_vals[:, 0] 
    y = curvature_vals[:, 1]  
    k = curvature_vals[:, 2]  

    dx = np.gradient(x)
    dy = np.gradient(y)
    dkdx = np.gradient(k, dx)  
    dkdy = np.gradient(k, dy) 
    dk = np.gradient(k)

    ds = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
    ds_mid = np.concatenate(([ds[0]], (ds[:-1] + ds[1:]) / 2, [ds[-1]]))
    dk_ds = dk / ds_mid

    print(f"dk size: {dk_ds.shape}")
    print(f"curv vals shape: {curvature_vals.shape}")
    print(dk**2)
    
    return dk, dk**2 #dkdx, dkdy
