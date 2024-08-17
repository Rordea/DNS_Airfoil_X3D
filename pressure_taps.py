# pressure_taps.py

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
import os

# polynomial_utils.py

import numpy as np
from sklearn.metrics import mean_squared_error



def press_taps(filename='NACA0018_upper.txt', degrees=range(5, 6), offset_distance=1.0, num_offset_points=10,initial_point_x=1,initial_point_y=1,initial_point_z=0):
    data = np.loadtxt(filename)
    x = data[:, 0]+initial_point_x
    y = data[:, 1]+initial_point_y

    best_degree = None
    min_error = float('inf')
    best_coefficients = None
    degree_errors = []

    for degree in degrees:
        coefficients = np.polyfit(x, y, degree)
        polynomial = np.poly1d(coefficients)
        y_fit = polynomial(x)
        mse = mean_squared_error(y, y_fit)
        degree_errors.append((degree, mse))
        
        if mse < min_error:
            min_error = mse
            best_degree = degree
            best_coefficients = coefficients

    best_polynomial = np.poly1d(best_coefficients)
    derivative = np.polyder(best_polynomial)

    offset_points = []
    for xi in np.linspace(min(x), max(x), num_offset_points):
        yi = best_polynomial(xi)
        slope = derivative(xi)
        normal_vector = np.array([-slope, 1])
        normal_vector /= np.linalg.norm(normal_vector)
        offset_point = np.array([xi, yi]) + offset_distance * normal_vector
        
        offset_points.append((offset_point[0], offset_point[1], initial_point_z) )

    offset_points = np.array(offset_points)

    results = {
        'best_degree': best_degree,
        'min_error': min_error,
        'best_coefficients': best_coefficients,
        'best_polynomial': best_polynomial,
        'degree_errors': degree_errors,
        'offset_points': offset_points
    }
    
    
    # print(f"Best Degree: {results['best_degree']}")
    # print(f"Minimum Mean Squared Error: {results['min_error']:.2f}")
    # print(f"The best polynomial coefficients are: {results['best_coefficients']}")
    # print(f"The best polynomial equation is: {results['best_polynomial']}\n")

    # # Plotting example (optional)
    # # for i, point in enumerate(results['offset_points']):
    # #     print(f"Pressure tap {i+1} coordinates: ({point[0]:.2f}, {point[1]:.2f}, {point[2]:.2f})")
    # #     probes.append( (point[0], point[1], point[2]) )
    # #probes.append(press_taps_coords['offset_points'])

    # x_fit = np.linspace(min(results['offset_points'][:, 0]), max(results['offset_points'][:, 0]), 100)
    # y_fit = results['best_polynomial'](x_fit)

    # plt.scatter(results['offset_points'][:, 0], results['offset_points'][:, 1], color='green', marker="+" ,label='Offset Points')
    # plt.plot(x_fit, y_fit, label=f'Best Degree {results["best_degree"]}, MSE: {results["min_error"]:.2f}', color='blue')
    # plt.legend()
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.title('Offset Points and Best Polynomial Fit')
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.show()  
    

    return results




# press_taps(filename='NACA0018_upper.txt', degrees=range(3,4), offset_distance=0.1, num_offset_points=10,initial_point_x=1,initial_point_y=1,initial_point_z=0)