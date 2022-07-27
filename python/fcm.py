class fcm:   
    """
    This is a firn compaction model and a few simple methods
    for performing parameter searches and plotting results, packaged as
    a class called fcm. 
    """

    def param_search_1D(self, param_to_vary, param_values, save_all_timesteps = False, **kwargs):
        '''
        This is intended as a convenient way of performing a parameter search. 
        It loops over a given parameter and runs the model for each value, then
        concatenates the results into one xarray dataset, with the parameter defined 
        as a new dimension and corresonding coordinate. 

        Inputs 
        param_to_vary -- a string of the name of the paraneter to perform the search over
        e.g., 'beta'

        param_values -- a list of numpy array of parameter values. These should be monotonically increasing.

        save_all_timesteps -- a logical determining if all the timesteps are to be saved. 
        If False (the default) then only the last time step of each simulation is saved. 

        **kwargs -- can be a series of keyword-value pairs overwriting the default values of 
        other model parameters, i.e. not the one whichthis 1D search is over, but any of the 
        other whose default values can be seen at the start of the 'run' function below. 

        Returns: a dataset in self.results which contains a new dimension for the parameter to 
        vary and has all the results concatenated along that dimension.

        Usage:
        simplest example:
        >>> import fcm 
        >>> sim = fcm.fcm()
        >>> sim.param_search_1D('beta', (0.5, 1))  
        The line above uses default values for all values other than beta, for which it uses 
        values 0.5 and 1. 
        
        >>> import fcm 
        >>> sim2 = fcm.fcm()
        >>> sim2.param_search_1D('r_s_dim', (2.5e-07, 3e-07, 4e-07), beta = 0.5)  
        The line above varies the dimensional grain size at the surface, and uses a non-default value for beta.

        '''
        
        import xarray as xr
        import numpy as np
        result_datasets = []                            # create an empty list to contain all the separate results xarray.datasets
        for x in param_values:                          # loop over the values of the parameter 
            print('Running model with ' + str(param_to_vary) + ' = ' + str(x))
            keywords_to_pass = {param_to_vary: x}       # construct a dictionary to pass to the solver with the appropriate parameter name and value
            self.run(print_messages = False, **keywords_to_pass, **kwargs)      # run the solver
            if save_all_timesteps:                  
                result_datasets.append(self.results)    # iteratively append the dataset into the list created above, if save_all_timesteps = True, then keep all the timesteps
            else:
                result_datasets.append(self.results.isel(t=-1))                     # iteratively append the dataset into the list created above, if save_all_timesteps = False, then just keep the last time step
        new_coords = xr.DataArray(np.array(param_values), dims=(param_to_vary))     # create an xarrat.dataarray contianing the values of the parameter for use in the line below when making a new dimension and coordinate in the results xarray. 
        self.results = xr.concat(result_datasets, new_coords)                       # concatenate all the datsets along a new dimension named after the parameter we are varying. Also add the values in which were put into the xarray in the line above.  

    def run_single(self, **kwargs):  
        self.run(**kwargs)



    def run(self, **kwargs):
        '''
        
        Usage: 
        >>> import importlib
        >>> sim = fcm.fcm()
        >>> sim.run()
        
        View the xarray.dataset of results as follows:
        >>> test.results


        '''
        
        import numpy as np
        self.p = {}

        # 1. define default paramaters   
        self.p['sim_label'] = 'default_label'
        self.p['b0_mpy'] = 0.1
        self.p['beta'] = 1
        self.p['nu'] = 1
        self.p['T_s_dim'] = 253.15
        self.p['plotting'] = 1
        self.p['saving_xt'] = 1
        self.p['dz'] = 0.01
        self.p['save'] = 1
        self.p['sim_T'] = True
        self.p['sim_r'] = True
        self.p['PauseGrainEvolution_t'] = np.nan
        self.p['z0'] = 100
        self.p['r_s_dim'] = 2.5000e-07   # grain size at the surface (0.5 mm)**2 
        self.p['phi_s'] = 0.5   
        self.p['n'] = 1
        self.p['m'] = 1
        self.p['simDuration'] = 4
        self.p['scaleDuration'] = False
        self.p['interp_on_reg_z'] = False 
        self.p['print_messages']   = True




        # 2. Replace any parameters that are defined by the user (in kwargs) 
        for key, values in kwargs.items():
            if 'print_messages' in locals() and print_messages == True:
                print(' ---- overwriting ', key, ' with a user-defined value of', str(values), '.')
            self.p[key] = values     


        if self.p['print_messages']:
            print('*** Starting setup.')
 

        # 3. Define (or extract from self.p) the dimensional parameters of the system.
        b0_mpy = self.p['b0_mpy']        # ice equivalent accumulation rate [m / yr]
        T_s_dim = self.p['T_s_dim']      # upper surface temperature [K]
        z_0 = self.p['z0']               # initial column height [m]
        dz_dim = self.p['dz']*z_0        # dimensional numerical grid spacing [m]
        r2_s_dim = self.p['r_s_dim']     # upper surface grain size [m**2] (0.5 mm)**2 
        r2_f = 0.01**2                        # maximum grain size [m**2] (1 cm)**2 
        phi_s = self.p['phi_s']          # upper surface porosity
        c_i = 2009                            # heat capacity [J / (kg K)]
        E_c = 60e3                            # compaction activation energy [J]
        E_g = 42e3                            # grain growth activation energy [J]
        k_c = 9.2e-9                          # flow parameter [kg?1 m3 s] from Arthern kc
        k_g = 1.3e-7/r2_f                     # grain growth rate constant [m**2 / s]
        m = self.p['m']                                  # porosity exponent
        n = self.p['n']                       # stress exponent
        kappa_0 = 2.1                         # thermal conductivity [W / (m K)], # arthern and wingham 1998 pg. 18
        G = 0.05                              # geothermal heat flux  [W/m**2]

        ## 4. Constants, scales and parameters
        ### 4.1 constants
        g = 9.81                           # acceleration due to gravity [m s**-2]
        spy = 24*365*3600                  # seconds per year
        R = 8.31                           # ideal gas constant
        rho_i = 918
        self.p['rho_i'] = rho_i                        # ice density [kg / m**3]

        ### 4.2 Scaling parameters.
        b_0 = b0_mpy/spy 
        h_0 = z_0 
        r2_0 = (h_0*k_g*r2_f/b_0)*np.exp(-E_g/(R*T_s_dim)) 
        t_0 = h_0/b_0 
        T_0 = G*z_0/kappa_0 
        sigma_0 = g*rho_i*h_0 
        w_0 = b_0                 # scale of the vertical velocity is the accumulation rate

        # save scales for use in plotting later
        self.p['b_0'] = b_0
        self.p['r2_0'] = r2_0
        self.p['t_0'] = t_0
        self.p['sigma_0'] = sigma_0
        self.p['w_0'] = w_0
        self.p['h_0'] = h_0
        self.p['T_0'] = T_0

        ### 4.3 Non-dimensional parameters.
        lambda_c = E_c*T_0/(R*T_s_dim**2.0) 
        lambda_g = E_g*T_0/(R*T_s_dim**2.0) 
        gamma = (sigma_0/4.0)**(1-n)    # factor to account for non-linear rheology, under the default value of n = 1, gamma = 1 and it has no effect on the value of Ar
        Ar = r2_0/(k_c*t_0*sigma_0*np.exp(-E_c/(R*T_s_dim)))/gamma  # Arthern number
        Fl = h_0*G/(kappa_0*T_0) 
        Pe = rho_i*c_i*b_0*h_0/kappa_0 
        beta = self.p['beta'] 
        delta = r2_0/r2_f 
        nu = self.p['nu']       # accumulation multiplier

        ## 5. Compute vertical velocity upper boundary conditon
        w_s = nu*beta/(1 - phi_s)  

        ## 6. Set up real space grid in height coordinates.
        z_init = np.arange(z_0,-dz_dim,-dz_dim)
        N = z_init.size

        ### Normalized depth coordinates.
        z_h = np.flip(z_init)/z_0 

        ## 6.1. time span
        simDuration = self.p['simDuration']
        beta = self.p['beta']

        if self.p['scaleDuration']:
            self.p['t_span'] = np.array([0, simDuration/beta])
        else:
            self.p['t_span']  = np.array([0, simDuration])


        ## 7. Initial conditions
        Ly0 = len(z_h)*4+1

        y0 = np.zeros(Ly0,dtype=np.float64)
        phi_init = y0[:-1:4]
        r2_hat_init = y0[1:-1:4]
        T_hat_init = y0[2:-1:4]
        A_hat_init = y0[3:-1:4]
        H_init = y0[Ly0-1:Ly0] 

        ### Initial porosity
        phi_init[:] = (1-z_h)*phi_s 

        ### Dimensional grain size squared.
        if self.p['sim_r']:
            r2_hat_init[:] = r2_s_dim/r2_0 + z_h 
        else:
            r2_hat_init[:] = r2_s_dim/r2_0 + 0*z_h 

        ### Dimensionless temperature (mostly not used)
        T_hat_init[:] = np.zeros(N) 

        ### Dimensionless firn age.
        A_hat_init[:] = z_h 

        ### domain height inital condition     
        H_init[:] = z_0/h_0

        ## 8. Define gradient operator
        ### Finite difference gradient operator using two-point upwind scheme.
        D1 = self.upwind_difference_matrix(z_h[0], z_h[-1], N, 1) 

        # 9. Add the newly created values to self.p
        self.p['N'] = N 
        self.p['D1'] = D1
        self.p['spy'] = spy 
        self.p['lambda_c'] = lambda_c 
        self.p['lambda_g'] = lambda_g 
        self.p['Ar'] = Ar 
        self.p['delta'] = delta 
        self.p['PecletNumber'] = Pe 
        self.p['FluxNumber'] = Fl 
        self.p['ArthenNumber'] = Ar 
        self.p['y0'] = y0
        self.p['z_h'] = z_h
        self.p['w_s'] = w_s

        if self.p['print_messages']:
            print('*** Setup complete.')     

    
        from scipy import integrate
        from scipy import interpolate
        import numpy.matlib
        import numpy as np
        import time
        import warnings

        ### main integration step
        st = time.time()   # record start time
        if self.p['print_messages']:
            print('*** Starting integration.')
        sol = integrate.solve_ivp(self.eqns, self.p['t_span'], self.p['y0']);    # need to change back to self.eqns
        et = time.time()    # record end time
        elapsed_time = et - st    # get exectuion time
        if self.p['print_messages']:
            if sol.status == 0:
                print(f"*** Succesfully finished integration in  {elapsed_time:.3} seconds.")
            else:
                warnings.warn('*** Integrator failed to find a solution.')
        
        ### Collect the variables.
        phi = sol.y[:-1:4]
        r2 = sol.y[1:-1:4]
        T = sol.y[2:-1:4]
        A = sol.y[3:-1:4]
        Height = sol.y[-1]
        t = sol.t
        
        ### post-processing
        if self.p['print_messages']:
            print('*** Post-processing.')

        ### compute density from porosity
        rho = self.p['rho_i']*(1 - phi)        
                    
        ### Un-normalized depth coordinates
        Hg, zg = np.meshgrid(Height,self.p['z_h'])
        z = Hg*zg     # non-scaled depth (still nondimensional) eqn B1 from TC manuscript

        n = self.p['n']
        m = self.p['m']
        Ar = self.p['Ar']
        z_h = self.p['z_h']
        lambda_c = self.p['lambda_c']
        w_s = self.p['w_s']
        z0 = self.p['z0']

        ### compute velocity, mass, firn air content and firn thickenss
        W = np.empty_like(phi)         # velocity
        M = np.empty_like(Height)      # total mass
        z830 = np.empty_like(Height)   # firn thickness (depth to 830 kg/m^3)
        FAC = np.empty_like(Height)    # firn air content
        for i in range(len(t)):
        ### Compute stress    
            S_int = Height[i]*(1 - phi[:,i])
            Sigma = integrate.cumulative_trapezoid(S_int,z_h,initial=0)
        ### Compute velocity    
            W_int = -(Height[i]/Ar)*Sigma**n*phi[:,i]*m*np.exp(lambda_c*T[:,i])/r2[:,i]
            W[:,i] = integrate.cumulative_trapezoid(W_int, z_h, initial=0) + w_s
        ### Compute total mass in the column.
            M[i] = integrate.trapezoid(1 - phi[:,i], z[:,i]*z0)
        ### Compute firn air content
            FAC[i] = integrate.trapezoid(phi[:,i], z[:,i]*z0)
        ### Compute the firn thickness 
            f = interpolate.interp1d(phi[:,i],z[:,i], bounds_error=False)
            z830[i] = f(1-830/self.p["rho_i"])

        ### Create output xarray
        import xarray as xr

        out = xr.Dataset(
            data_vars=dict(
                phi=(["z_h", "t"], phi),
                r2=(["z_h", "t"], r2),
                rho=(["z_h", "t"], rho),
                A=(["z_h", "t"], A),
                T=(["z_h", "t"], T),
                w=(["z_h", "t"], W),
                h=(["t"], Height),
                M=(["t"], M),
                FAC=(["t"], FAC),
                z830=(["t"], z830),
            ),
            coords=dict(
                z_h=(["z_h"], self.p['z_h']),
                z=(["z_h", "t"], z),
                t=t,
            ),
            attrs=dict(simulation_parameters=self.p),
        )


        # add attributes to variables
        out.phi.attrs = dict(name="porosity", 
            long_name="nondimensional porosity",
            scale=1,
            scale_units="none")

        out.rho.attrs = dict(name="density", 
            long_name="density",
            units="m^3/s")

        out.r2.attrs = dict(name="grainsize", 
            long_name="nondimensional squared grain radius",
            scale=self.p["r2_0"],
            scale_units="m^2")

        out.A.attrs = dict(name="age", 
            long_name="nondimensional age",
            scale=self.p["t_0"],
            scale_units="s")

        out.T.attrs = dict(name="temperature", 
            long_name="nondimensional temperature",
            scale="to be added")

        out.h.attrs = dict(name="height", 
            long_name="nondimensional domain height",
            scale=self.p["z0"],
            scale_units="m")

        out.w.attrs = dict(name="velocity", 
            long_name="nondimensional velocity",
            scale=self.p['w_0'],
            scale_units="m/s")

        out.M.attrs = dict(name="mass", 
            long_name="total ice-equivelent mass",
            scale_units="m ice-equivelent",
            notes="This is the total mass in the domain expressed as meters ice-equivelent. It is computed as the depth integral of the 1-phi at each time step.")

        out.FAC.attrs = dict(name="firn air content", 
            long_name="firn air content",
            units="m",
            notes="Computed as the depth integral of the porosity. Note that if the firn is not completed compacted by the bottom of the domain (i.e. phi(z_h=1,t) > 0) then this FAC value does not accurately represent the total firn air content, as some air nominally exists beneath the bottom of the domain.")

        out.z830.attrs = dict(name="firn thickness", 
            long_name="nondimensional firn thickenss",
            scale=self.p['z0'],
            scale_units="m")

        # add attributes to coordinates
        out.z_h.attrs = dict(name="non_dim_depth",
            long_name="scaled and nondimensional depth",
            scale="H(t) * z_0",
            scale_units="m",
            notes="the scaling of this variable changes with time through the simulation. So z_h = 1 means a different dimensional depth at each time step. Specifically, z_h = 1 corresponds to H(t)*z_0. H(t) is one of the other variables, and z_0 is the dimensional scale for depth [m].")
            
        out.t.attrs = dict(name="time", 
            long_name="nondimensional time",
            scale=self.p["t_0"],
            scale_units="s")

        out.z.attrs = dict(name="depth", 
            long_name="nondimensional depth, irregular grid",
            scale=self.p['z0'],
            scale_units="m")

        self.results = out

        # interpolate onto regular vertical grid, z
        if self.p['interp_on_reg_z'] == True:
            if self.p['print_messages']:
                print('**** interpolating onto regular vertical grids')
            self.interp_regular_z(var_name = 'phi')
            self.interp_regular_z(var_name = 'rho')
            self.interp_regular_z(var_name = 'r2')
            self.interp_regular_z(var_name = 'A')
            self.interp_regular_z(var_name = 'T')
            self.interp_regular_z(var_name = 'w')

    # plot depth profiles
    def profiles(self,var_to_plot='phi', time_to_plot=1000):
        '''
        Simple plots of the depth variation of model variables. 

        var_to_plot is the variable you want to plot. 
        This has to be one which varies with depth: 
        e.g., phi, r2, rho, A, or T. The default is the porosity, phi.

        time_to_plot is the time steps to plot, in nondimensional time. 
        This can be a single value or a list. The default is 100. Because this 
        is much larger than the duration of the simulation, this means that the 
        default is the final time step. 

        For example, if the simulation have been run as follows:
        >>> import fcm
        >>> sim = fcm.fcm()
        >>> sim.integrate()

        You can plot a single profile of porosity with
        >>> sim.profiles('phi',0.4)
        
        or you can plot several profiles from different time slices with
        >>> sim.profiles('phi',[0, 0.1, 4])
        or 
        >>> sim.profiles('A',[0, 0.1, 4])

        If no time steps correspon exactly to the requests simulation times, 
        the nearest time step will be selected. 

        '''
        
        da = self.results[var_to_plot]
        da.sel(t=time_to_plot, method='nearest').plot.line(x='z_h') 

# plot z-t diagrams
    def z_t_plot(self,var_to_plot='phi'):
        '''
        Depth-time plots.  

        var_to_plot is the variable you want to plot. 
        This has to be one which varies with depth and time: 
        e.g., phi, r2, rho, A, or T. The default is the porosity, phi.


        For example, if the simulation have been run as follows:
        >>> import fcm
        >>> sim = fcm.fcm()
        >>> sim.integrate()

        you can plot a depth-time plot with
        >>> sim.z_t_plot('phi')

        If no time steps correspond exactly to the requests simulation times, 
        the nearest time step will be selected. 

        '''
        
        da = self.results[var_to_plot]
        da.plot(yincrease = False) # profiles of last time slice

    def time_series(self,var_to_plot='H'):
        '''
        Time series of variables that do not vary with depth: H, FAC, z830, M  

        Inputs: 
        var_to_plot is the variable you want to plot. 
        This has to be one which does not vary with depth: 
        e.g., phi, r2, rho, A, or T. The default is the porosity, phi.


        Returns: 
        Creates a plot of the a time series of the prescribed variable.

        Usage:
        For example, if the simulation have been run as follows:
        >>> import fcm
        >>> sim = fcm.fcm()
        >>> sim.integrate()

        you can plot a time series plot with
        >>> sim.time_series('FAC')
        '''          

        da = self.results[var_to_plot]
        da.plot.line(x='t')


    def interp_at_depth(self, depth_to_interp = 20, var_name = 'phi'):
        '''
        Interpolating variables at a prescribed depth for all time steps. 

        Inputs:
        depth_to_interp -- the depth in meters you want the variable interpolated at (a number)
        var_name -- the name of the variable you want interplolated (str)

        Returns:
        Adds a variable to the results xarray dataset with a name name from the two inputs, for example, results.phi_20

        Usage:
        If the simulation has been run as follows:
        >>> import fcm
        >>> sim = fcm.fcm()
        >>> sim.integrate()

        Use
        >>> sim.interp_at_depth(depth_to_interp = 30, var_name = 'A')
        '''
        
        import warnings
        warnings.warn('This method is outdated. You are better off applying fcm.interp_regular_z to the variable you are interested in, to interpolate it on to a regular vertical grid.')

        from scipy import interpolate
        import numpy as np
        import xarray as xr

        name_for_new_var = var_name + '_' + str(depth_to_interp)

        interpolated_values = np.empty_like(self.results.t.values)
        for i in range(len(self.results.t.values)):
            f = interpolate.interp1d(self.results.z.isel(t=i).values*self.p['z0'],self.results[var_name].isel(t=i).values, bounds_error = False)
            interpolated_values[i] =f(depth_to_interp)


        self.results[name_for_new_var] = self.results.h.copy(data=interpolated_values)



    def interp_regular_z(self, var_name = 'phi'):
        '''
        Interpolate variable on to a regular - i.e. a grid which is same at all time steps. 
        
        Inputs:
        var_name -- the name of the variable you want interplolated (str)

        Returns:
        Adds a variable to the results xarray dataset with '_r' affixed to the end  (standing for 'regular')
        of the variable name, e.g., phi_r.

        It also adds a new dimension-coordinate, z_r, which is the nondimensional depth on a regular grid. 
        This is the same for every time step so the coordinate is only a function of z. 
        
        Usage:
        If the simulation has been run as follows:
        >>> import fcm
        >>> sim = fcm.fcm()
        >>> sim.integrate()

        Use
        >>> sim.interp_regular_z(var_name = 'A')
        This will produce a new variable called A_r and a new dimension cooridnate called z_r. 
        '''
      
        from scipy import interpolate
        import numpy as np
        import xarray as xr
       
        z_q = np.linspace(0,1.2,round(self.p['N']*1.2))   # query points for interpolation
        Nz = len(z_q)
        Nt = len(self.results.t.values)
        interpolated_values = np.empty((Nz,Nt))

        for i in range(Nt):
            
            #z_h = self.results.z_h.values                 # scaled depth
            #h =   self.results.h.isel(t=i).values         # domain length
            #z =   z_h*h                                   # nonscaled depth (still nondimensional): eqn B1 from TC manuscript
            z =   self.results.z.isel(t=i).values 
            Y =   self.results[var_name].isel(t=i).values # variable on old grid

            f = interpolate.interp1d(z,Y,bounds_error = False)
            interpolated_values[:,i] =f(z_q)

        
        name_for_new_var = var_name + '_r'

        t_attrs = self.results.t.attrs    # save the attributes of t, because the line below deletes them. 

        self.results[name_for_new_var] = xr.DataArray(
            data=interpolated_values,
            dims=["z_r", "t"],
            coords=dict(
                z_r=z_q))


        self.results.t.attrs = t_attrs

        new_name_for_xarray = self.results[var_name].attrs['name'] + ' (regular)'
        new_long_name_for_xarray = self.results[var_name].attrs['long_name'] + ' (regular grid)'

        self.results[name_for_new_var].attrs = dict(name=new_name_for_xarray, 
            long_name=new_long_name_for_xarray,
            scale=self.p['z0'],
            scale_units="m")
        
        
        self.results.z_r.attrs = dict(name="depth", 
            long_name="nondimensional depth, regular grid",
            scale=self.p['h_0'],
            scale_units="m")

    # Model equations      
    def eqns(self,t, y):
        from scipy import integrate
        import numpy as np

        Ly = len(y)
        
        N = self.p['N']   
        z_h = self.p['z_h']
        D1 = self.p['D1']
        lambda_g  = self.p['lambda_g']
        lambda_c  = self.p['lambda_c']
        delta = self.p['delta']
        Ar = self.p['ArthenNumber']
        w_s = self.p['w_s']
        beta = self.p['beta']
        n = self.p['n']
        m = self.p['m']

        ### Collect the simulation variables.
        phi = y[:-1:4]
        r2 = y[1:-1:4]
        T = y[2:-1:4]
        A = y[3:-1:4]
        H = y[-1]

        ### create dydt vector
        dydt = np.empty_like(y)

        ### create views of dydt
        dphidt = dydt[:-1:4]
        dr2dt = dydt[1:-1:4]
        dTdt = dydt[2:-1:4]
        dAdt = dydt[3:-1:4]
        dHdt = dydt[Ly-1:Ly]

        ### Compute gradients.
        dphidz = np.matmul(D1, phi) 
        dr2dz  = np.matmul(D1, r2)  
        dAdz   = np.matmul(D1, A)  

        ### Compute stress.
        s_int = H*(1 - phi) 
        sigma = integrate.cumulative_trapezoid(s_int,z_h, initial=0) 

        ### Compute the ice velocity.
        v_int = -(H/Ar)*sigma**n*phi**m*np.exp(lambda_c*T)/r2 
        w = integrate.cumulative_trapezoid(v_int, z_h, initial=0) + w_s 

        ### Column height.
        dHdt[:] = w[-1] - beta/(1 - phi[-1])    

        ### Change in porosity.
        dphidt[:] = (1/H)*(D1@((1 - phi)*w) + dHdt*z_h*dphidz) 

        ### Change in square of the grain size.
        if self.p['sim_r']:   # (only if sim_r ==1) 
            dr2dt[:] = (1/H)*(dHdt*z_h - w)*dr2dz + (1 - delta*r2)*np.exp(lambda_g*T) 
        else:
            dr2dt[:] = 0*r2 

        ### Change in temperature.
        dTdt[:] = 0 

        ### Change in age.
        dAdt[:] = 1 + (1/H)*(dHdt*z_h - w)*dAdz 

        ### Upper surface boundary conditions.
        dphidt[0] = 0 
        dr2dt[0] = 0 
        dAdt[0] = 0 

                
        return dydt

    # upwind_difference_matrix
    def upwind_difference_matrix(self,z0,zL,n,v):
        # Adapted from Kerschbaum, Simon. (2020). Backstepping Control of Coupled Parabolic Systems with Varying Parameters: A Matlab Library (1.0). Zenodo. https://doi.org/10.5281/zenodo.4274740       
        
        
        from scipy.sparse import spdiags
        import numpy as np
        dz=(zL-z0)/(n-1)

        # (1) finite difference approximation for positive v     
        if v > 0:
            # sparse discretization matrix      
            # interior points      
            e = np.ones(n)
            D = spdiags(np.stack([-e, e]), np.array([-1, 0]), n, n).toarray()
            # boundary point      
            D[0,0:2] = np.array([-1, +1])

        # (2) finite difference approximation for negative v
        if v < 0:
            # interior points      
            e = np.ones(n)
            D = spdiags(np.stack([-e, e]), np.array([0, 1]), n, n).toarray()
            # boundary point      
            D[n-1,(n-2):(n)] = np.array([-1, +1])                        
        
        return D/dz