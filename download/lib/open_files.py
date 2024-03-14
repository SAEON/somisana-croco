from contextlib import contextmanager


@contextmanager
def open_files(*file_paths, **kwargs):
    file_handles = []
    try:
        for file_path in file_paths:
            file_handle = open(file_path, **kwargs)
            file_handles.append(file_handle)
        yield file_handles
    finally:
        for file_handle in file_handles:
            file_handle.close()
