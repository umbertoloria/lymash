def ksliding(source, window_size, string=True):
	offset = 0
	windows = []
	while offset + window_size <= len(source):
		window = source[offset:offset + window_size]
		if not string:
			if 0 not in window:
				windows.append(window)
		else:
			windows.append(window)
		offset += 1
	return windows
