import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load data from CSV
data = pd.read_csv('results.csv', names=['Implementation', 'Matrix Size', 'No of Threads', 'Time (ms)'])

# Convert 'Time (ms)' to seconds
data['Time (s)'] = data['Time (ms)'] / 1000

# Extract sequential times for speedup calculation
sequential_times = data[data['Implementation'] == 'SEQUEN'].set_index('Matrix Size')['Time (s)']

# Calculate speedup
def calculate_speedup(row):
    if row['Implementation'] != 'SEQUEN':
        seq_time = sequential_times.get(row['Matrix Size'], None)
        if seq_time is not None:
            return seq_time / row['Time (s)']
    return None

# Apply speedup calculation to parallel implementations
data['Speedup'] = data.apply(calculate_speedup, axis=1)

# Function to plot speedup bar chart
def plot_speedup(implementation_data, title):
    matrix_sizes = sorted(implementation_data['Matrix Size'].unique())
    thread_counts = sorted(implementation_data['No of Threads'].unique())
    
    # Bar width and positions
    bar_width = 0.1
    index = np.arange(len(matrix_sizes))
    
    # Plot bars
    plt.figure(figsize=(15, 7))
    for i, thread_count in enumerate(thread_counts):
        speedups = [implementation_data[(implementation_data['Matrix Size'] == size) & (implementation_data['No of Threads'] == thread_count)]['Speedup'].values[0] for size in matrix_sizes]
        plt.bar(index + i * bar_width, speedups, bar_width, label=f'{thread_count} Threads')

    plt.xlabel('Matrix Size (n)')
    plt.ylabel('Speedup')
    plt.title(title)
    plt.xticks(index + bar_width / 2, matrix_sizes)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Filter data for OpenMP and PThread and plot the speedup bar chart
openmp_data = data[data['Implementation'] == 'OPENMP']
pthread_data = data[data['Implementation'] == 'PTHREAD']

plot_speedup(openmp_data, 'Speedup with OpenMP')
plot_speedup(pthread_data, 'Speedup with Pthreads')
