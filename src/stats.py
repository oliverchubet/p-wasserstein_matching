import pstats

# Load the stats file
stats = pstats.Stats('wass.prof')

# Clean up filenames (optional, makes it easier to read)
stats.strip_dirs()

# Sort by cumulative time (usually the most useful metric)
stats.sort_stats('cumulative')

# Print the top 20 lines
stats.print_stats(20)
