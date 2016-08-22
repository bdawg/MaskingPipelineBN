    # =========================================================================
    # =========================================================================

    def chi2_3d(self,smin,smax,anglemin,anglemax,cmin,cmax,
                   nsep=100,nangle=360,nc=100):

        ''' Generates a chi2 grid in all three variables. Cannibalised from
        Frantz' code to do this in pyapm/plot_search. Goes over three variables,
        so takes longer than fit_map, which is just for looking at projections
        quickly.

        The indices are [i,j,k] over [sep,angle,contrast].

        NOTE FOR ANTHONY: This is all object oriented, so it's a method in a
        class that holds all of the data. Remove all references to "self" to
        make this take ordinary variables. I assume you have a phase_binary
        function somewhere...'''

        filter = self.info['filter']

        chi2 = np.zeros((nsep, nangle, nc))

        # define our grid points

        seps = smin + (smax-smin) * np.arange(nsep)/nsep
        angles  = anglemin + (anglemax-anglemin) * np.arange(nangle)/nangle
        cons = cmin + (cmax-cmin) * np.arange(nc)/nc

        # Calculate null hypothesis chi-squared

        nullchi2 = np.sum((self.kp_signal/self.kp_error)**2)

        for i,sep in enumerate(seps):
            print i, sep
            for j,angle in enumerate(angles):
                for k,c in enumerate(cons):
                    test = phase_binary(self.uv[:,0], self.uv[:,1],
                                        filter, [sep,angle+90.0-self.info['orient'],c])
                    modl_ker = np.dot(self.KerPhi, self.RED*test)
                    #vari[i,j,k] = np.std(modl_ker-self.kp_signal)
                    chi2[i,j,k] = np.sum(((modl_ker-self.kp_signal)/self.kp_error)**2)

        #wherevar = np.where(vari == vari.min())

        wh = np.where(chi2 == chi2.min())

        bestsep, bestangle, bestcon = seps[wh[0][0]], angles[wh[1][0]], cons[wh[2][0]]

        #print 'Best-fit separation = %.2f, position angle = %.2f, contrast = %.2f' % (bestsep,bestangle,bestcon)

        #print 'Best chi2 obtained = %.3f' % chi2.min()

        

        return chi2
