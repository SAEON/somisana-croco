import psutil
import os

def get_memory_usage():
    process = psutil.Process(os.getpid()) # Get process details
    mem_info = process.memory_info() # Get process memory info
    return bytes_to_human(mem_info.rss) # returns the resident set size in bytes

# Scale bytes to a readable formats
def bytes_to_human(bytes, suffix="B"):
    factor = 1024
    for unit in ["", "K", "M", "G", "T", "P"]:
        if bytes < factor:
            return f"{bytes:.2f}{unit}{suffix}"
        bytes /= factor
