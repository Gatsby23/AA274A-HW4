import numpy as np

EPSILON_OMEGA = 1e-3

def compute_dynamics(xvec, u, dt, compute_jacobians=True):
    """
    Compute Turtlebot dynamics (unicycle model).

    Inputs:
                     xvec: np.array[3,] - Turtlebot state (x, y, theta).
                        u: np.array[2,] - Turtlebot controls (V, omega).
        compute_jacobians: bool         - compute Jacobians Gx, Gu if true.
    Outputs:
         g: np.array[3,]  - New state after applying u for dt seconds.
        Gx: np.array[3,3] - Jacobian of g with respect to xvec.
        Gu: np.array[3,2] - Jacobian of g with respect to u.
    """
    ########## Code starts here ##########
    # TODO: Compute g, Gx, Gu
    # HINT: To compute the new state g, you will need to integrate the dynamics of x, y, theta
    # HINT: Since theta is changing with time, try integrating x, y wrt d(theta) instead of dt by introducing om
    # HINT: When abs(om) < EPSILON_OMEGA, assume that the theta stays approximately constant ONLY for calculating the next x, y
    #       New theta should not be equal to theta. Jacobian with respect to om is not 0.
    theta = xvec[2]
    x = xvec[0]
    y = xvec[1]

    V = u[0]
    w = u[1]

    s_w = w
    inv_w = None

    g = None
    Gx = None
    Gu = None

    if abs(w) < EPSILON_OMEGA:
        #TODO: Add constant-angle calculations here
        sin_t = np.sin(theta) + np.sin(theta + w*dt)
        cos_t = np.cos(theta) + np.cos(theta + w*dt)

        g_lst = [x + V*0.5*cos_t*dt, y + V*0.5*sin_t*dt, theta + w*dt]
        g = np.array(g_lst)

        Gx_lst = [[1, 0, -0.5*V*sin_t*dt], [0, 1, 0.5*V*cos_t*dt], [0, 0, 1]]
        Gx = np.array(Gx_lst)


        Gu_lst = [[0.5*cos_t*dt, -0.5*V*np.sin(theta + w*dt)*dt*dt], [0.5*sin_t*dt, 0.5*V*np.cos(theta + w*dt)*dt*dt], [0, dt]]
        Gu = np.array(Gu_lst)
    else:
        inv_w = 1.0 / w
        n_theta = theta + s_w * dt
        upper_sin = np.sin(n_theta)
        lower_sin = np.sin(theta)

        upper_cos = np.cos(n_theta)
        lower_cos = np.cos(theta)


        j_theta = theta + w * dt
        j_upper_sin = np.sin(j_theta)
        j_lower_sin = np.sin(theta)
        j_upper_cos = np.cos(j_theta)
        j_lower_cos = np.cos(theta)

        #NOTE: if w < EPSILON_OMEGA, we need to calculate everything totally differently.  We can't integrate over theta now.  It is just V*cos(theta)*dt. This also has different jacobians
        g_lst = [x + V*inv_w*(upper_sin - lower_sin), y + V*inv_w*(-upper_cos + lower_cos), theta + w*dt]
        g = np.array(g_lst)

        Gx_lst = [[1, 0, V*inv_w*(j_upper_cos - j_lower_cos)], [0, 1, V*inv_w*(j_upper_sin - j_lower_sin)], [0, 0, 1]]

    
        x_dw = V*(-1.0 * (inv_w**2) * (j_upper_sin - j_lower_sin) + inv_w * (j_upper_cos * dt))
        y_dw = V*(-1.0 * (inv_w**2) * (-j_upper_cos + j_lower_cos) + inv_w * (j_upper_sin * dt))
        Gu_lst = [[inv_w*(j_upper_sin - j_lower_sin), x_dw], [inv_w*(-j_upper_cos + j_lower_cos), y_dw], [0, dt]]

        Gu = np.array(Gu_lst)
        Gx = np.array(Gx_lst)
    ########## Code ends here ##########

    if not compute_jacobians:
        return g

    #print ("{0}, {1}, {2}".format(g, Gx, Gu))
    return g, Gx, Gu
'''def compute_dynamics(xvec, u, dt, compute_jacobians=True):
    theta = xvec[2]
    x = xvec[0]
    y = xvec[1]

    V = u[0]
    w = u[1]
    s_w = w
    if abs(w) < EPSILON_OMEGA:
        s_w = 0

    n_theta = theta + w * dt
    upper_sin
    
    g_lst = []
    if not compute_jacobians:
        return g'''
    
    
    


def transform_line_to_scanner_frame(line, x, tf_base_to_camera, compute_jacobian=True):
    """
    Given a single map line in the world frame, outputs the line parameters
    in the scanner frame so it can be associated with the lines extracted
    from the scanner measurements.

    Input:
                     line: np.array[2,] - map line (alpha, r) in world frame.
                        x: np.array[3,] - pose of base (x, y, theta) in world frame.
        tf_base_to_camera: np.array[3,] - pose of camera (x, y, theta) in base frame.
         compute_jacobian: bool         - compute Jacobian Hx if true.
    Outputs:
         h: np.array[2,]  - line parameters in the scanner (camera) frame.
        Hx: np.array[2,3] - Jacobian of h with respect to x.
    """
    alpha, r = line

    ########## Code starts here ##########
    x_cam, y_cam, th_cam = tf_base_to_camera
    x_base, y_base, th_base = x

    tf_robot_to_world = np.array([[np.cos(th_base), -np.sin(th_base), x_base],
                                  [np.sin(th_base),  np.cos(th_base), y_base], 
                                  [              0,                0,      1]])
    
    x_cam_world, y_cam_world, th_cam_world = tf_robot_to_world.dot(np.array([x_cam, y_cam, 1]))
    
    h = np.array([alpha - th_base - th_cam, 
                  r - x_cam_world*np.cos(alpha) - y_cam_world*np.sin(alpha)])

    Hx_23 =  (y_cam*np.cos(alpha)-x_cam*np.sin(alpha))*np.cos(th_base)
    Hx_23 += (x_cam*np.cos(alpha)+y_cam*np.sin(alpha))*np.sin(th_base)

    Hx = np.array([[             0,              0,    -1], 
                   [-np.cos(alpha), -np.sin(alpha), Hx_23]])
    ########## Code ends here ##########

    if not compute_jacobian:
        return h

    return h, Hx


def normalize_line_parameters(h, Hx=None):
    """
    Ensures that r is positive and alpha is in the range [-pi, pi].

    Inputs:
         h: np.array[2,]  - line parameters (alpha, r).
        Hx: np.array[2,n] - Jacobian of line parameters with respect to x.
    Outputs:
         h: np.array[2,]  - normalized parameters.
        Hx: np.array[2,n] - Jacobian of normalized line parameters. Edited in place.
    """
    alpha, r = h
    if r < 0:
        alpha += np.pi
        r *= -1
        if Hx is not None:
            Hx[1,:] *= -1
    alpha = (alpha + np.pi) % (2*np.pi) - np.pi
    h = np.array([alpha, r])

    if Hx is not None:
        return h, Hx
    return h
