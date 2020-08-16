def ksliding(source, window_size):
	offset = 0
	windows = []
	while offset + window_size <= len(source):
		window = source[offset:offset + window_size]
		windows.append(window)
		offset += 1
	return windows
