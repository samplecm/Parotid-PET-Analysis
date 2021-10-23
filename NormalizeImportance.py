import numpy as np 

def SubImportance():
    subImportance = np.array([-0.0005628706, 0.0074248731, -0.0001937888, 0.027761328, 0.0049947, 0.00796205, 0.00191083, 0.0165397])
    normalizedImportance = (subImportance - np.amin(subImportance)) / (np.amax(subImportance) - np.amin(subImportance))
    return normalizedImportance.tolist()

