import random

def get_random_color(pastel_factor=0.5):
	return [(x+pastel_factor)/(1.0+pastel_factor) for x in [random.uniform(0,1.0) for i in [1,2,3]]]

def color_distance(c1,c2):
	return sum([abs(x[0]-x[1]) for x in zip(c1,c2)])

def generate_new_color(existing_colors, pastel_factor=0.5):
	max_distance = None
	best_color = None
	for i in range(0, 100):
		color = get_random_color(pastel_factor = pastel_factor)
		if not existing_colors:
			return color
		best_distance = min([color_distance(color,c) for c in existing_colors])
		if not max_distance or best_distance > max_distance:
			max_distance = best_distance
			best_color = color
	return best_color

def generate_N_pastels(N=10, color_format='hex', pastel_factor=0.5):
	colors = []
	for i in range(0,N):
		colors.append(generate_new_color(colors, pastel_factor))
	if color_format == 'hex':
		colors = [('#%02x%02x%02x' % tuple(map(lambda k: 255*k, x))).upper() for x in colors]
	elif color_format == 'dec':
		colors = [tuple(map(lambda k: 255*k, x)) for x in colors]
	return colors

if __name__ == '__main__':

  #To make your color choice reproducible, uncomment the following line:
  #random.seed(10)
	outfile = open("random_colors_1000.txt", "wb")
	colors =  generate_N_pastels(1000)
	outfile.write("\n".join(colors))
	outfile.close()
