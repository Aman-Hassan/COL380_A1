import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import re

# Matrix size
matrix_size = input("Enter matrix size: ")

# Number of threads
threads = [1, 2, 4, 8, 16]

# Dataframe to store results
df = pd.DataFrame(columns=['Threads', 'LUpthreads Efficiency', 'LUomp Efficiency'])

# Function to run C program
def run_program(program, threads, matrix_size):
    # Run the program with different numbers of threads
    times = []
    for t in threads:
        print(f"Running {program} with {t} threads and matrix size {matrix_size}...")
        result = subprocess.run(["./bin/" + program, str(matrix_size), str(t)], capture_output=True, text=True)
        # Extract execution time from output
        time = float(re.search(r"Execution time: (\d+\.\d+) seconds", result.stdout).group(1))
        print(f"Execution time: {time} seconds")
        times.append(time)

    # Calculate efficiencies
    efficiencies = [times[0] / (t * time) for t, time in zip(threads, times)]
    return efficiencies

print("Starting the process...")

subprocess.run(["make", "all"])
# Run LUpthreads and LUomp
df['Threads'] = threads
print("Running LUpthreads...")
df['LUpthreads Efficiency'] = run_program("LUpthreads", threads, matrix_size)
print("Running LUomp...")
df['LUomp Efficiency'] = run_program("LUomp", threads, matrix_size)

print("Generating the graphs...")

# Plot LUpthreads graph
plt.figure(figsize=(8, 6))
plt.plot(df['Threads'], df['LUpthreads Efficiency'], label='LUpthreads')
plt.xlabel('Number of threads')
plt.ylabel('Parallel Efficiency')
plt.legend()
plt.savefig('./plots/LUpthreads.png')
plt.clf()

# Plot LUomp graph
plt.figure(figsize=(8, 6))
plt.plot(df['Threads'], df['LUomp Efficiency'], label='LUomp')
plt.xlabel('Number of threads')
plt.ylabel('Parallel Efficiency')
plt.legend()
plt.savefig('./plots/LUomp.png')
plt.clf()

# Plot both together
plt.figure(figsize=(8, 6))
plt.plot(df['Threads'], df['LUpthreads Efficiency'], label='LUpthreads')
plt.plot(df['Threads'], df['LUomp Efficiency'], label='LUomp')
plt.xlabel('Number of threads')
plt.ylabel('Parallel Efficiency')
plt.legend()
plt.savefig('./plots/Both.png')
plt.clf()

subprocess.run(["make", "clean"])
# Display table
print("Here is the table of results:")
print(df)

print("Process completed.")