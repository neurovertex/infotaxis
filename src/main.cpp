#include "infotaxis.hpp"
#include <iostream>
#include <algorithm>
#include <list>
#include <cmath>
#include <png++/png.hpp>

using namespace std;
using namespace png;
using namespace infotaxis;

rgb_pixel get_colour(double absv, double relv) {
	double value = min((int)(relv*32)*8, 255);
	int b = (int)(value * cos(absv*M_PI/2)),
			r = (int)(value * sin(absv*M_PI/2)), 
			g = b;
	//cout << "relv : " << relv << ", absv : " << absv << ", rgb : "<< r <<":"<< g <<":"<< b << endl;
	return rgb_pixel(r, g, b);
}

typedef struct simultrace {
	int x, y, detections;
} Trace;

struct simul {
	int curx, cury;
	double time;
	InfotaxisGrid *grid;
	list<Trace> backtrace;
};

void write_image(string filename, double *array, int w, int h, struct simul *simul, const int ratio) {
	const int imgw = w*ratio, imgh = h*ratio;
	image<rgb_pixel> image(imgw, imgh);
	double mi = INFINITY, ma = -INFINITY;
	for (int i = 0; i < h; i ++)
		for (int j = 0; j < w; j ++) {
			double raw = array[i*w + j], val = (raw);
			if (std::isinf(raw))
				cerr << "Infinite : "<< j << ":"<< i << endl;
			else if (std::isnan(raw))
				cerr << "NaN : "<< j << ":"<< i << endl;
			//cout << raw << endl;
			if (!std::isnan(val) && !std::isinf(val)) {
				mi = min(mi, val);
				ma = max(ma, val);
			} else
				val = min(max(val, mi), ma);
			array[i*w + j] = val;
		}

	double scale = (ma != mi) ? ma-mi : 1;

	for (int i = 0; i < imgh; i ++)
		for (int j = 0; j < imgw; j ++) {
	//		cout << j << ":" << i <<" =>\t" << array[i*w/ratio + j/ratio] << endl;
			image[imgh-1-i][j] = get_colour(array[i/ratio*w + j/ratio], (array[i/ratio*w + j/ratio]-mi)/scale);
		}
	cout << filename << " : max : "<< ma <<", min : "<< mi << endl;
	if (simul != NULL) {
		for (Trace t : simul->backtrace) {
			rgb_pixel color = (t.detections == 0) ? rgb_pixel(64, 128, 128) : rgb_pixel(255 - (128/t.detections), 64, 64);
			for (int j = 1; j < max(ratio-1, 2); j ++)
				for (int i = 1; i < max(ratio-1, 2); i ++) {
					image[imgh-1-(t.y * ratio + j)][t.x * ratio + i] = color;
				}
		}

		int x = simul->curx * ratio + ratio/2;
		int y = simul->cury * ratio + ratio/2;
		for (int i = 0; i <= 10; i ++) {
			if (imgh-y+4-i >= 0 && imgh-y+4-i < imgh) {
				if (x-5+i >= 0 && x-5+i < imgw)
					image[imgh-y+4-i][x-5+i] = rgb_pixel(255, 128, 0);
				
				if (x+5-i >= 0 && x+5-i < imgw)
					image[imgh-y+4-i][x+5-i] = rgb_pixel(255, 128, 0);
			}
		}
	}

	image.write(filename);
}

void write_image(string filename, double *array, int w, int h) {
	write_image(filename, array, w, h, NULL, 1);
}

int main() {
	const int w=100, h=100, x0 = w/2, y0 = h/2+10, ttl = 400;
	const double diff=1, rate=1, windagl = -M_PI/2, windvel = 2, a = .1, dt = 1, resolution = 0.05;
	InfotaxisGrid grid(w, h, diff, rate, windvel, windagl, ttl, a, resolution);
	struct simul simul = {w/4, h/4, 0, &grid};
	char *filename = new char[64];
	int cnt = 0;
  	std::default_random_engine gen;
	gen.seed(time(NULL));
	/*grid.updateProbas(simul.curx++, simul.cury, 0, simul.time++);
	write_image("output0.png", grid[0], w, h);*/
	cout << "Proba sous source : "<< grid.concentration(x0, y0-10, x0, y0) << endl;

	while ((simul.curx != x0 || simul.cury != y0) && cnt < 1000) {
		double mean = grid.encounterRate(simul.curx, simul.cury, x0, y0) * dt;
		cout << cnt <<" : "<< simul.curx <<":"<< simul.cury <<" => mean "<< mean;
		poisson_distribution<int> pdist(mean);
		gen(); // Creating a variate_generator with gen seeds it from gen's value without modifying gen.
		//Thus we need to iterate gen to avoid reseeding the generator with the same value every iteration
		int detects = pdist(gen);
		grid.updateProbas(simul.curx, simul.cury, detects, simul.time);
		cout <<", detections : "<< detects << ", entropy : "<< grid.entropy() << endl;
		simul.time += dt;
		Direction optimal = grid.getOptimalMove(simul.curx, simul.cury, dt);
		//cout << ", optimal move : "<< DIRECTIONS[optimal+1] << "("<< optimal <<")" << endl;
		sprintf(filename, "pictures/iteration%04d.png", cnt++);
		write_image(filename, grid[0], w, h, &simul, 5);
		simul.backtrace.push_front({simul.curx, simul.cury, detects});
		go_to(optimal, simul.curx, simul.cury);
	}
	sprintf(filename, "pictures/iteration%04d.png", cnt);
	write_image(filename, grid[0], w, h, &simul, 5);
	if (simul.curx == x0 && simul.cury == y0)
		cout << "Source found after "<< cnt << " iterations at " << x0 <<":"<< y0 <<" !" << endl;
	else
		cout << "Source not found after "<< cnt <<" iterations"<< endl;
}


