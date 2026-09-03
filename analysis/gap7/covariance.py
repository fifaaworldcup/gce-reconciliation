import numpy as np
from pathlib import Path

def read_results(file_path):
  data = np.loadtxt(file_path)
  return data

def compute_covariance(data):
  covariance_matrix = np.cov(data, rowvar=False)
  return covariance_matrix

def main():
  file_path = Path("results.txt")

data = read_results(file_path)
covariance_matrix = compute_covariance(data)

print(covariance_matrix)

if _name_ == "_main_":
  main()

