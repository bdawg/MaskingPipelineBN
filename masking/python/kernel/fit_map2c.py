def fit_map2c(a,b,sep=None,angle=None,c1=None, c2 = None,
                sepmin = 40.0,sepmax = 1000.0,anglemin=0.0,
                anglemax=359.9, cmin = 1.01,cmax = 10.0,
                display= True,nclevs = 20,clevsmax=12,filename=None):
        ''' A version of fit_map that can handle two datasets in different
        colours to generate combined chi-squared maps. We enforce a single
        angular separation and position angle and allow c1 and c2 to vary.
        There are then '''

        filter = self.info['filter']

        # first make the parameter space

        nx,ny = 120,120


        #calculate a regular grid
        if sep is not None:
            key = 1
            ttl = 'Reduced chi-squared map, sep = ' + str(sep)
            xaxis = np.linspace(anglemin,anglemax,num=nx)
            xlbl = 'Position Angle (degrees)'
            xmin, xmax = anglemin, anglemax
            yaxis = np.linspace(cmin,cmax,num=ny)
            ylbl = 'Contrast'
            ymin, ymax = cmin, cmax
            xy = np.meshgrid(xaxis,yaxis)
        elif angle is not None:
            key = 2
            ttl = 'Reduced chi2 map, position angle = '+str(angle)
            angle = angle
            xaxis = np.linspace(sepmin,sepmax,num=nx)
            xlbl = 'Separation (mas)'
            xmin, xmax = sepmin,sepmax
            yaxis = np.linspace(cmin,cmax,num=ny)
            ylbl = 'Contrast'
            ymin, ymax = cmin, cmax
        elif contrast is not None:
            key = 3
            ttl = 'Reduced chi2 map, contrast = ' + str(contrast)
            xaxis = np.linspace(anglemin,anglemax,num = nx)
            xlbl = 'Position Angle (degrees)'
            xmin, xmax = anglemin, anglemax
            yaxis = np.linspace(sepmin,sepmax,num=ny)
            ylbl = 'Separation (mas)'
            ymin, ymax = sepmin, sepmax
        else:
            print "fit_map2c passed invalid argument(s)"

        # generate a uv phase map for each of the parameter sets

        xx,yy = np.meshgrid(xaxis,yaxis)

        #calculate models

        if key ==1:
            models = [[phase_binary(self.uv[:,0],self.uv[:,1],filter,[sep,xx[j][i]+90.0 -self.info['orient'],yy[j][i]])
                       for i in range(0,nx)] for j in range(0,ny)]
        elif key ==2:
            models = [[phase_binary(self.uv[:,0],self.uv[:,1],filter,[xx[j][i],angle+90.0 -self.info['orient'],yy[j][i]])
                       for i in range(0,nx)] for j in range(0,ny)]
        elif key ==3:
            models = [[phase_binary(self.uv[:,0],self.uv[:,1],filter,[yy[j][i],xx[j][i]+90.0 -self.info['orient'],contrast])
                       for i in range(0,nx)] for j in range(0,ny)]
            
        #calculate kerphases
            
        kphasemodels = [[np.dot(self.KerPhi,self.RED*models[j][i]) for i in range (0,nx)]
                        for j in range(0,ny)]

        #calculate departures from the model and the null hypothesis

        diff = [[np.divide(np.subtract(self.kp_signal,kphasemodels[j][i])**2,
                           self.kp_error**2) for i in range(0,nx)] for j in range(0,ny)]

        null = np.divide(self.kp_signal**2,self.kp_error**2)
        
        #calculate chi2 for the model and the null hypothesis

        chi2 = [[np.sum(diff[j][i]) for i in range(0,nx)] for j in range(0,ny)]
        
        chi2null= np.sum(null)

        # calculate reduced chi2

        dof = float(self.nkphi - 4)

        chi2 = np.array(chi2)/dof

        chi2null = np.array(chi2null)/dof

        # calculate significance

        bestfit = chi2.min()

        significance = np.sqrt(chi2null - bestfit)

        # find the best-fit coordinates

        bestindices = np.unravel_index(chi2.argmin(),chi2.shape)

        bestx, besty = xx[bestindices], yy[bestindices]

        if display is True: #default, plot it
            print 'Best chi-squared = ', bestfit

            print 'Null hypothesis chi-squared = ', chi2null

            print 'Significance estimated at ', significance

            print 'Best fit coordinates = ', bestx, besty

            #create plot

            clevs = chi2.min()*np.linspace(1,clevsmax,nclevs)

            plt.contourf(xx,yy,chi2,levels = clevs,cmap=plt.cm.bone,
                         antialiased = True)
            plt.axis('tight')
            plt.xlabel(xlbl)
            plt.ylabel(ylbl)
            plt.title(ttl)
            plt.colorbar()
            
            if filename is not None:
                plt.savefig(filename+'.png',format='png')
            else:
                plt.show()
            
            return None
        else: # return parameters in a string p = [bestfit,significance,bestx,besty]
            p = [bestfit,significance,bestx,besty]
            return p
