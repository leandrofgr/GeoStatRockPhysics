import scipy.io as sio
import numpy as np

teste = sio.loadmat('Sw24_inverted.mat')
data = teste['Sw24_inverted']
matlab_data = np.array(data)
np.save("Sw24_inverted.npy", np.ascontiguousarray(matlab_data))