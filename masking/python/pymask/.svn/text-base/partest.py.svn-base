import numpy as np
import playdoh
import matplotlib.pyplot as plt

def wastetime(seed):
	rand = np.random.randn(500)
	return np.std(rand)

seeds = np.linspace(0,1,100)

#now do in parallel

stuff = playdoh.map(wastetime,seeds)