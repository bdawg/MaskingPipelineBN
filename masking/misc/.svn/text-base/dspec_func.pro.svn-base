;This is for fitting to dark speckle residuals

FUNCTION dspec_func, r, theta, p, dparms
      
               RETURN, p[0] + p[1]*cos(theta) + p[2]*sin(theta) +p[3]*cos(2*theta) + p[4]*sin(2*theta)+$
                 (p[5] + p[6]*cos(theta) + p[7]*sin(theta) + p[8]*cos(2*theta) + p[9]*sin(2*theta))/(r+1.0) +$
                 (p[10] + p[11]*cos(theta) + p[12]*sin(theta) + p[13]*cos(2*theta) + p[14]*sin(2*theta))*r
END
