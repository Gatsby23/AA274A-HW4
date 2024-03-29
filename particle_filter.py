import numpy as np
import scipy.linalg  # You may find scipy.linalg.block_diag useful
import scipy.stats  # You may find scipy.stats.multivariate_normal.pdf useful
import turtlebot_model as tb

EPSILON_OMEGA = 1e-3

class ParticleFilter(object):
    """
    Base class for Monte Carlo localization and FastSLAM.

    Usage:
        pf = ParticleFilter(x0, R)
        while True:
            pf.transition_update(u, dt)
            pf.measurement_update(z, Q)
            localized_state = pf.x
    """

    def __init__(self, x0, R):
        """
        ParticleFilter constructor.

        Inputs:
            x0: np.array[M,3] - initial particle states.
             R: np.array[2,2] - control noise covariance (corresponding to dt = 1 second).
        """
        self.M = x0.shape[0]  # Number of particles
        self.xs = x0  # Particle set [M x 3]
        self.ws = np.repeat(1. / self.M, self.M)  # Particle weights (initialize to uniform) [M]
        self.R = R  # Control noise covariance (corresponding to dt = 1 second) [2 x 2]

    @property
    def x(self):
        """
        Returns the particle with the maximum weight for visualization.

        Output:
            x: np.array[3,] - particle with the maximum weight.
        """
        idx = self.ws == self.ws.max()
        x = np.zeros(self.xs.shape[1:])
        x[:2] = self.xs[idx,:2].mean(axis=0)
        th = self.xs[idx,2]
        x[2] = np.arctan2(np.sin(th).mean(), np.cos(th).mean())
        return x

    def transition_update(self, u, dt):
        """
        Performs the transition update step by updating self.xs.

        Inputs:
            u: np.array[2,] - zero-order hold control input.
            dt: float        - duration of discrete time step.
        Output:
            None - internal belief state (self.xs) should be updated.
        """
        ########## Code starts here ##########
        # TODO: Update self.xs.
        # Hint: Call self.transition_model().
        # Hint: You may find np.random.multivariate_normal useful.
        dt = np.array(dt).astype(np.float64)
        eps = np.random.multivariate_normal(np.zeros(u.shape[0]), self.R*dt, size=self.M) # Gaussian input noise
        us = u[None,:] + eps
        self.xs = self.transition_model(us, dt)
        ########## Code ends here ##########

    def transition_model(self, us, dt):
        """
        Propagates exact (nonlinear) state dynamics.

        Inputs:
            us: np.array[M,2] - zero-order hold control input for each particle.
            dt: float         - duration of discrete time step.
        Output:
            g: np.array[M,3] - result of belief mean for each particle
                               propagated according to the system dynamics with
                               control u for dt seconds.
        """
        raise NotImplementedError("transition_model must be overridden by a subclass of EKF")

    def measurement_update(self, z_raw, Q_raw):
        """
        Updates belief state according to the given measurement.

        Inputs:
            z_raw: np.array[2,I]   - matrix of I columns containing (alpha, r)
                                     for each line extracted from the scanner
                                     data in the scanner frame.
            Q_raw: [np.array[2,2]] - list of I covariance matrices corresponding
                                     to each (alpha, r) column of z_raw.
        Output:
            None - internal belief state (self.x, self.ws) is updated in self.resample().
        """
        raise NotImplementedError("measurement_update must be overridden by a subclass of EKF")

    def resample(self, xs, ws):
        """
        Resamples the particles according to the updated particle weights.

        Inputs:
            xs: np.array[M,3] - matrix of particle states.
            ws: np.array[M,]  - particle weights.

        Output:
            None - internal belief state (self.xs, self.ws) should be updated.
        """
        r = np.random.rand() / self.M

        ########## Code starts here ##########
        # TODO: Update self.xs, self.ws.
        # Note: Assign the weights in self.ws to the corresponding weights in ws
        #       when resampling xs instead of resetting them to a uniform
        #       distribution. This allows us to keep track of the most likely
        #       particle and use it to visualize the robot's pose with self.x.
        # Hint: To maximize speed, try to implement the resampling algorithm
        #       without for loops. You may find np.linspace(), np.cumsum(), and
        #       np.searchsorted() useful. This results in a ~10x speedup.
        ws_cumsum = np.cumsum(ws)
        ws_total = ws_cumsum[-1] # last element in ws_cumsum is the overall total
        m = np.linspace(0, self.M, self.M, False) # 0, 1, 2, ... M-1
        u = ws_total * (r + m/self.M)
        particle_idxs = np.searchsorted(ws_cumsum, u)
        self.xs = xs[particle_idxs]
        self.ws = ws[particle_idxs]
        ########## Code ends here ##########

    def measurement_model(self, z_raw, Q_raw):
        """
        Converts raw measurements into the relevant Gaussian form (e.g., a
        dimensionality reduction).

        Inputs:
            z_raw: np.array[2,I]   - I lines extracted from scanner data in
                                     columns representing (alpha, r) in the scanner frame.
            Q_raw: [np.array[2,2]] - list of I covariance matrices corresponding
                                     to each (alpha, r) column of z_raw.
        Outputs:
            z: np.array[2I,]   - joint measurement mean.
            Q: np.array[2I,2I] - joint measurement covariance.
        """
        raise NotImplementedError("measurement_model must be overridden by a subclass of EKF")


class MonteCarloLocalization(ParticleFilter):

    def __init__(self, x0, R, map_lines, tf_base_to_camera, g):
        """
        MonteCarloLocalization constructor.

        Inputs:
                       x0: np.array[M,3] - initial particle states.
                        R: np.array[2,2] - control noise covariance (corresponding to dt = 1 second).
                map_lines: np.array[2,J] - J map lines in columns representing (alpha, r).
        tf_base_to_camera: np.array[3,]  - (x, y, theta) transform from the
                                           robot base to camera frame.
                        g: float         - validation gate.
        """
        self.map_lines = map_lines  # Matrix of J map lines with (alpha, r) as columns
        self.tf_base_to_camera = tf_base_to_camera  # (x, y, theta) transform
        self.g = g  # Validation gate
        super(self.__class__, self).__init__(x0, R)

    def transition_model(self, us, dt):
        """
        Unicycle model dynamics.

        Inputs:
            us: np.array[M,2] - zero-order hold control input for each particle.
            dt: float         - duration of discrete time step.
        Output:
            g: np.array[M,3] - result of belief mean for each particle
                               propagated according to the system dynamics with
                               control u for dt seconds.
        """

        ########## Code starts here ##########
        # TODO: Compute g.
        # Hint: We don't need Jacobians for particle filtering.
        # Hint: A simple solution can be using a for loop for each partical
        #       and a call to tb.compute_dynamics
        # Hint: To maximize speed, try to compute the dynamics without looping
        #       over the particles. If you do this, you should implement
        #       vectorized versions of the dynamics computations directly here
        #       (instead of modifying turtlebot_model). This results in a
        #       ~10x speedup.
        # Hint: This faster/better solution does not use loop and does 
        #       not call tb.compute_dynamics. You need to compute the idxs
        #       where abs(om) > EPSILON_OMEGA and the other idxs, then do separate 
        #       updates for them
        g = np.zeros((self.M, 3))
        theta = self.xs[...,2]
        x = self.xs[...,0]
        y = self.xs[...,1]
        V = us[...,0]
        w = us[...,1]
        s_w = w

        sin_t = np.sin(theta) + np.sin(theta + w*dt)
        cos_t = np.cos(theta) + np.cos(theta + w*dt)
        g1 = np.stack([x + V*0.5*cos_t*dt, y + V*0.5*sin_t*dt, theta + w*dt], -1)

        inv_w = 1.0 / np.maximum(np.abs(w),EPSILON_OMEGA)*np.sign(w)
        n_theta = theta + s_w * dt
        upper_sin = np.sin(n_theta)
        lower_sin = np.sin(theta)
        upper_cos = np.cos(n_theta)
        lower_cos = np.cos(theta)
        g2 = np.stack([x + V*inv_w*(upper_sin - lower_sin), y + V*inv_w*(-upper_cos + lower_cos), theta + w*dt], -1)

        g = np.where(np.abs(w[:,None])<EPSILON_OMEGA, g1, g2)

        # for i in range(self.M):
        #     x = self.xs[i]
        #     u = us[i]
        #     h = tb.compute_dynamics(x, u, dt, compute_jacobians=False)
        #     g[i] = h
        ########## Code ends here ##########

        return g

    def measurement_update(self, z_raw, Q_raw):
        """
        Updates belief state according to the given measurement.

        Inputs:
            z_raw: np.array[2,I]   - matrix of I columns containing (alpha, r)
                                     for each line extracted from the scanner
                                     data in the scanner frame.
            Q_raw: [np.array[2,2]] - list of I covariance matrices corresponding
                                     to each (alpha, r) column of z_raw.
        Output:
            None - internal belief state (self.x, self.ws) is updated in self.resample().
        """
        xs = np.copy(self.xs)
        ws = np.zeros_like(self.ws)

        ########## Code starts here ##########
        # TODO: Compute new particles (xs, ws) with updated measurement weights.
        # Hint: To maximize speed, implement this without looping over the
        #       particles. You may find scipy.stats.multivariate_normal.pdf()
        #       useful.
        # Hint: You'll need to call self.measurement_model()
        vs, Q = self.measurement_model(z_raw, Q_raw)
        ws = scipy.stats.multivariate_normal.pdf(vs, mean=None, cov=Q)
        ########## Code ends here ##########

        self.resample(xs, ws)

    def measurement_model(self, z_raw, Q_raw):
        """
        Assemble one joint measurement and covariance from the individual values
        corresponding to each matched line feature for each particle.

        Inputs:
            z_raw: np.array[2,I]   - I lines extracted from scanner data in
                                     columns representing (alpha, r) in the scanner frame.
            Q_raw: [np.array[2,2]] - list of I covariance matrices corresponding
                                     to each (alpha, r) column of z_raw.
        Outputs:
            z: np.array[M,2I]  - joint measurement mean for M particles.
            Q: np.array[2I,2I] - joint measurement covariance.
        """
        vs = self.compute_innovations(z_raw, np.array(Q_raw))

        ########## Code starts here ##########
        # TODO: Compute Q.
        # Hint: You might find scipy.linalg.block_diag() useful
        Q = scipy.linalg.block_diag(*Q_raw)
        ########## Code ends here ##########

        return vs, Q

    def compute_innovations(self, z_raw, Q_raw):
        """
        Given lines extracted from the scanner data, tries to associate each one
        to the closest map entry measured by Mahalanobis distance.

        Inputs:
            z_raw: np.array[2,I]   - I lines extracted from scanner data in
                                     columns representing (alpha, r) in the scanner frame.
            Q_raw: np.array[I,2,2] - I covariance matrices corresponding
                                     to each (alpha, r) column of z_raw.
        Outputs:
            vs: np.array[M,2I] - M innovation vectors of size 2I
                                 (predicted map measurement - scanner measurement).
        """
        def angle_diff(a, b):
            a = a % (2. * np.pi)
            b = b % (2. * np.pi)
            diff = a - b
            if np.size(diff) == 1:
                if np.abs(a - b) > np.pi:
                    sign = 2. * (diff < 0.) - 1.
                    diff += sign * 2. * np.pi
            else:
                idx = np.abs(diff) > np.pi
                sign = 2. * (diff[idx] < 0.) - 1.
                diff[idx] += sign * 2. * np.pi
            return diff

        ########## Code starts here ##########
        # TODO: Compute vs (with shape [M x I x 2]).
        # Hint: Simple solutions: Using for loop, for each particle, for each 
        #       observed line, find the most likely map entry (the entry with 
        #       least Mahalanobis distance).
        # Hint: To maximize speed, try to eliminate all for loops, or at least
        #       for loops over J. It is possible to solve multiple systems with
        #       np.linalg.solve() and swap arbitrary axes with np.transpose().
        #       Eliminating loops over J results in a ~10x speedup.
        #       Eliminating loops over I results in a ~2x speedup.
        #       Eliminating loops over M results in a ~5x speedup.
        #       Overall, that's 100x!
        # Hint: For the faster solution, you might find np.expand_dims(), 
        #       np.linalg.solve(), np.meshgrid() useful.

        J = self.map_lines.shape[1]
        I = z_raw.shape[1]

        # hs = self.compute_predicted_measurements()
        # vs = np.empty((self.M, z_raw.shape[1], 2))
        # for i, h in enumerate(hs):
        #     for j in range(I):
        #         z = z_raw[:, j]
        #         innovations = np.zeros_like(h.T)
        #         innovations[:, 0] = angle_diff(z[0], h[0])
        #         innovations[:, 1] = z[1] - h[1]
        #         # malhanobis distances
        #         Q = Q_raw[j,:]
        #         Q_inv = np.linalg.inv(Q)
        #         distances = np.array([np.matmul(v.T, np.matmul(Q_inv, v)) for v in innovations])
        #         index = np.argmin(distances)
        #         vs[i, j] = innovations[index]

        hs = self.compute_predicted_measurements().transpose(0, 2, 1) # (M, J, 2)
 
        z_matrix = z_raw.T[None, None, :, :] # (1, 1, I, 2)
        h_matrix = hs[:, :, None, :] # (M, J, 1, 2)
        
        # Vectorized angle_diff
        z_alpha = z_matrix[..., 0] % (2.*np.pi) # (M, J, I)
        h_alpha = h_matrix[..., 0] % (2.*np.pi)
        angle_diffs = z_alpha - h_alpha 
        idxs = np.pi < angle_diffs
        signs = 2 * (angle_diffs[idxs] < 0) - 1
        angle_diffs[idxs] += signs * 2. * np.pi
        innovations_alpha = angle_diffs

        innovations_r = z_matrix[..., 1] - h_matrix[..., 1]
        innovations = np.stack((innovations_alpha, innovations_r), axis=3)

        v = innovations[..., None] # (M, J, I, 2, 1)
        Q_inv = np.linalg.inv(Q_raw)[None, None, :, :, :] # (1, 1, I, 2, 2)

        d_matrix = np.matmul(np.matmul(v.transpose(0, 1, 2, 4, 3), Q_inv), v) # (M, J, I, 1, 1)
        d_matrix = d_matrix.reshape((self.M, J, I))  # (M, J, I)

        d_argmin = np.argmin(d_matrix, axis=1)[:, None, :, None] # (M, 1, I, 1)
        vs = np.take_along_axis(innovations, d_argmin, axis=1) # (M, 1, I, 2)
        vs = vs.reshape((self.M, I, 2)) # (M, I, 2)
        ########## Code ends here ##########

        # Reshape [M x I x 2] array to [M x 2I]
        return vs.reshape((self.M,-1))  # [M x 2I]

    def compute_predicted_measurements(self):
        """
        Given a single map line in the world frame, outputs the line parameters
        in the scanner frame so it can be associated with the lines extracted
        from the scanner measurements.

        Input:
            None
        Output:
            hs: np.array[M,2,J] - J line parameters in the scanner (camera) frame for M particles.
        """
        ########## Code starts here ##########
        # TODO: Compute hs.
        # Hint: We don't need Jacobians for particle filtering.
        # Hint: Simple solutions: Using for loop, for each particle, for each 
        #       map line, transform to scanner frmae using tb.transform_line_to_scanner_frame()
        #       and tb.normalize_line_parameters()
        # Hint: To maximize speed, try to compute the predicted measurements
        #       without looping over the map lines. You can implement vectorized
        #       versions of turtlebot_model functions directly here. This
        #       results in a ~10x speedup.
        # Hint: For the faster solution, it does not call tb.transform_line_to_scanner_frame()
        #       or tb.normalize_line_parameters(), but reimplement these steps vectorized.       
        
        J = self.map_lines.shape[1]
        # hs = np.empty((self.M, 2, J))
        # for i, particle in enumerate(self.xs): # for each particle
        #     h = np.zeros_like(self.map_lines)
        #     for j, line in enumerate(self.map_lines.T): # for each map line
        #       h_col = tb.transform_line_to_scanner_frame(line, particle, self.tf_base_to_camera, False)
        #         h[:,j] = tb.normalize_line_parameters(h_col)
        #     hs[i] = h
        # return hs

        hs = self.map_lines.T # (J, 2)
        alpha, r = hs[:, 0], hs[:, 1]

        # Vectorized transform_line_to_scanner_frame
        x_cam, y_cam, th_cam = self.tf_base_to_camera
        x_base, y_base, th_base = self.xs.T

        tf_robot_to_world = np.array([[np.cos(th_base), -np.sin(th_base), x_base],
                                      [np.sin(th_base),  np.cos(th_base), y_base], 
                                      [              0,                0,      1]])
        
        x_cam_world, y_cam_world, th_cam_world = tf_robot_to_world.dot(np.array([x_cam, y_cam, 1]))


        alpha_cam = alpha[None, :] - th_base[:, None] - th_cam
        r_cam = (r[None, :] - x_cam_world[:, None]*np.cos(alpha)[None, :]
                            - y_cam_world[:, None]*np.sin(alpha)[None, :])
        
        # Vectorized normalize_line_parameters
        idxs = r_cam < 0
        alpha_cam[idxs] += np.pi
        r_cam[idxs] *= -1
        alpha_cam = (alpha_cam + np.pi) % (2*np.pi) - np.pi
        hs = np.array([alpha_cam, r_cam]).transpose(1, 0, 2) # (n, 2, n_lin)
        ########## Code ends here ##########

        return hs

