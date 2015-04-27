#include "infotaxis.hpp"
#include <algorithm>
#include <list>

using namespace std;
using namespace png;
using namespace infotaxis;

typedef struct simultrace {
	int x, y, detections;
} Trace;

struct simul {
	int curx, cury;
	double time;
	InfotaxisGrid *grid;
	list<Trace> backtrace;
};

void writeProbas(string filename, struct simul &simul, const int ratio, const int *x0 = NULL, const int *y0 = NULL) {
	const int imgw = simul.grid->getWidth() * ratio, imgh = simul.grid->getHeight() * ratio;
	image<rgb_pixel> image((size_t) imgw, (size_t) imgh);

	if (x0 == NULL || y0 == NULL)
		simul.grid->writeProbabilityField(image, ratio);
	else
		simul.grid->writeMeanStationaryField(image, *x0, *y0, ratio);

	for (Trace t : simul.backtrace) {
		rgb_pixel color = (t.detections == 0) ? rgb_pixel(64, 128, 128) : rgb_pixel(
				(byte) (255 - (128 / t.detections)), 64, 64);
		for (int j = 1; j < max(ratio - 1, 2); j++)
			for (int i = 1; i < max(ratio - 1, 2); i++) {
				image[imgh - 1 - (t.y * ratio + j)][t.x * ratio + i] = color;
			}
	}
	int x = simul.curx * ratio + ratio / 2;
	int y = simul.cury * ratio + ratio / 2;
	for (int i = 0; i <= 10; i++) {
		if (imgh - y + 4 - i >= 0 && imgh - y + 4 - i < imgh) {
			if (x - 5 + i >= 0 && x - 5 + i < imgw)
				image[imgh - y + 4 - i][x - 5 + i] = rgb_pixel(255, 128, 0);

			if (x + 5 - i >= 0 && x + 5 - i < imgw)
				image[imgh - y + 4 - i][x + 5 - i] = rgb_pixel(255, 128, 0);
		}
	}

	image.write(filename);
}

int main() {
	const int w = 100, h = 100, x0 = w / 2, y0 = h / 2 + 10, ttl = 400;
	const double diff = 1, rate = 1, windagl = -M_PI / 2, windvel = 2, a = .1, dt = 1, resolution = 0.05;
	InfotaxisGrid grid(w, h, diff, rate, windvel, windagl, ttl, a, resolution);
	struct simul simul = {w / 4, h / 4, 0, &grid};
	char *filename = new char[64];
	int cnt = 0;
	std::default_random_engine gen;
	gen.seed(time(NULL));

	while ((simul.curx != x0 || simul.cury != y0) && cnt < 1000) {
		double mean = grid.encounterRate(simul.curx, simul.cury, x0, y0) * dt;
		cout << cnt << " : " << simul.curx << ":" << simul.cury << " => mean " << mean;
		poisson_distribution<int> pdist(mean);
		gen(); // Creating a variate_generator with gen seeds it from gen's value without modifying gen.
		//Thus we need to iterate gen to avoid reseeding the generator with the same value every iteration
		int detects = pdist(gen);
		grid.updateProbas(simul.curx, simul.cury, detects, simul.time);
		cout << ", detections : " << detects << ", entropy : " << grid.entropy() << endl;
		simul.time += dt;
		Direction optimal = grid.getOptimalMove(simul.curx, simul.cury, dt);
		//cout << ", optimal move : "<< DIRECTIONS[optimal+1] << "("<< optimal <<")" << endl;
		sprintf(filename, "pictures/iteration%04d.png", cnt++);
		writeProbas(filename, simul, 5, NULL);
		simul.backtrace.push_front({simul.curx, simul.cury, detects});
		go_to(optimal, simul.curx, simul.cury);
	}
	sprintf(filename, "pictures/iteration%04d.png", cnt);
	writeProbas(filename, simul, 5, NULL, NULL);
	writeProbas("pictures/meanfield.png", simul, 5, &x0, &y0);
	if (simul.curx == x0 && simul.cury == y0)
		cout << "Source found after " << cnt << " iterations at " << x0 << ":" << y0 << " !" << endl;
	else
		cout << "Source not found after " << cnt << " iterations" << endl;
}
