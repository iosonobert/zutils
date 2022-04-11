import matplotlib.pyplot as plt
import numpy as np



def get_ellipse2D(var_maj, var_min, theta, n):
    """
    Get coordinates of an n point ellipse on an x-y plane with major axis at anb angle of theta. var_maj 
    is the length of the major axis and var_min is the length of the minor axis. 
    
    """
    t = np.linspace(0, 2*np.pi, n)
    x = var_maj*np.cos(t)*np.cos(theta) - var_min*np.sin(t)*np.sin(theta)
    y = var_maj*np.cos(t)*np.sin(theta) + var_min*np.sin(t)*np.cos(theta)

    return x, y

def PCA_2D(X, Y, plotopt=False):
    """
    PCA on a plane (2D)
    
    Originally from matlab code of Jeff Book's masterclass, modified by Tamara Schlosser, and then by Andrew Zulberti to 
    be more similar to the Emery and Thompson notation. 
    
    TO DO: generalise to 3D
    
    Inputs:
        X - timeseries 1 - on arbitrary axis x
        Y - timeseries 2 - on axis perpendicular to x
        
    Outputs:
        theta - radian angle of PC measured counterclockwise from axis x
        x     - x coordinates of ellipse
        y     - y coordinates of ellipse
    """
    
    # %% Now attempting the notation in Emery and Thompson

    X = X[~np.isnan(X)]
    Y = Y[~np.isnan(Y)]

    mean_x = np.mean(X)
    mean_y = np.mean(Y)

    u_1 = X - mean_x
    u_2 = Y - mean_y

    c = np.cov(u_1, u_2)
    cov_12 = c[0, 1]
    var_1  = c[0, 0]
    var_2  = c[1, 1]

    # % det|C - lambda*I| = 0
    # % using this in the quadratic equation gives
    lambda_1 = (1/2)*((var_1+var_2)+np.sqrt((var_1-var_2)**2 + 4*cov_12**2))
    lambda_2 = (1/2)*((var_1+var_2)-np.sqrt((var_1-var_2)**2 + 4*cov_12**2))

    theta_p = (1/2)*np.arctan2(2*cov_12, (var_1 - var_2))

    sqrt_lambda_1 = np.sqrt(lambda_1)
    sqrt_lambda_2 = np.sqrt(lambda_2)

    x, y = get_ellipse2D(2*sqrt_lambda_1, 2*sqrt_lambda_2, theta_p, 100)
#     x, y = get_ellipse2D(sqrt_lambda_1, sqrt_lambda_2, theta_p, 100)

    if plotopt:
        
        n = len(X)

        print('Ellipse axis: {} deg'.format(theta_p*180/np.pi))
        
        plt.figure()
        plt.plot(X, Y, '.')
        plt.plot(x, y, 'r-')
        plt.title('n points: {} | Ellipse axis: {:0.2f} deg'.format(n, theta_p*180/np.pi))
    
    return theta_p, x, y