import random
from scipy.stats import truncnorm
from tqdm import tqdm
from point import Point
import numpy as np

def generate_points(n, d, distribution):
    """
    Generates the point sets with the corresponding cost matrix
    :param compute_matrix: a function that computes the cost matrix
    """
    masses_a = [1 for _ in range(n)]
    masses_b = [1 for _ in range(n)]

    if distribution == "Normal_same":
        lower_bound = 0
        upper_bound = 1
        mean = 0.4
        std_dev = 0.1
        
        na = (lower_bound - mean) / std_dev
        nb = (upper_bound - mean) / std_dev
        
        A = [
            #Point(list(truncnorm.rvs(na, nb, loc=mean, scale=std_dev, size=d)), i, masses_a[i])
            Point(list(truncnorm.rvs(na, nb, loc=mean, scale=std_dev, size=d)), None)
            for i in range(n)
        ]
        B = [
            #Point(list(truncnorm.rvs(na, nb, loc=mean, scale=std_dev, size=d)), n + i, masses_b[i])
            Point(list(truncnorm.rvs(na, nb, loc=mean, scale=std_dev, size=d)), None)
            for i in range(n)
        ] 
    elif distribution == "Normal_different":
        # Define the parameters
        lower_bound = 0
        upper_bound = 1
        mean = 0.3
        std_dev = 0.3
        
        na = (lower_bound - mean) / std_dev
        nb = (upper_bound - mean) / std_dev
        
        A = [
            #Point(list(truncnorm.rvs(na, nb, loc=mean, scale=std_dev, size=d)), i, masses_a[i])
            Point(list(truncnorm.rvs(na, nb, loc=mean, scale=std_dev, size=d)), None)
            for i in range(n)
        ]
        mean = 0.7
        na = (lower_bound - mean) / std_dev
        nb = (upper_bound - mean) / std_dev
        B = [
            #Point(list(truncnorm.rvs(na, nb, loc=mean, scale=std_dev, size=d)), n + i, masses_b[i])
            Point(list(truncnorm.rvs(na, nb, loc=mean, scale=std_dev, size=d)), None)
            for i in range(n)
        ] 
    elif distribution == "Random_plane":
        while True:
            # Step 1: Randomly generate a basis vector for the plane.
            basis_1 = np.random.rand(d)
            basis_2 = np.random.rand(d)

            # Normalize the basis vectors.
            basis_1 /= np.linalg.norm(basis_1)
            basis_2 /= np.linalg.norm(basis_2)
            
            cosine_similarity = np.dot(basis_1, basis_2)  # Cosine of the angle
            if cosine_similarity <= 1/2:  # Ensure the angle >= 60 degrees
                break
        
        A = [
            #Point(list((random.random() / 2) * basis_1 + (random.random() / 2) * basis_2), i, masses_a[i])
            Point(list((random.random() / 2) * basis_1 + (random.random() / 2) * basis_2), None)
            for i in range(n)
        ]
        B = [
            #Point(list((random.random() / 2) * basis_1 + (random.random() / 2) * basis_2), n + i, masses_b[i])
            Point(list((random.random() / 2) * basis_1 + (random.random() / 2) * basis_2), None)
            for i in range(n)
        ]
    else:
        A = [
            #Point([float(random.random()) for _ in range(d)], i, masses_a[i])
            Point([float(random.random()) for _ in range(d)], None)
            for i in range(n)
        ]
        B = [
            #Point([float(random.random()) for _ in range(d)], n + i, masses_b[i])
            Point([float(random.random()) for _ in range(d)], None)
            for i in range(n)
        ]

    return A, B, masses_a, masses_b
