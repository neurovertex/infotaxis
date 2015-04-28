#include "infotaxis.hpp"
#include <algorithm>
#include <list>
#include <boost/program_options.hpp>

#ifdef INFOTAXIS_NO_PROFILING
#define PROFILE(cmd)
#define MSECONDS(time) (double)(.0)
#else
#include <time.h>
#define PROFILE(cmd) if (profiling) cmd;
#define MSECONDS(time) ((time * timeratio * 1000.))
#endif

#define MAX_ITERATIONS 2048

using namespace std;
using namespace png;
using namespace infotaxis;
namespace po = boost::program_options;

typedef struct simultrace {
	int x, y, detections;
} Trace;

struct simul {
	int curx, cury;
	double time;
	InfotaxisGrid *grid;
	list<Trace> backtrace;
};

void writeProbas(string filename, struct simul &simul, const int ratio, const int *x0 = NULL, const int *y0 = NULL);



int main(int argc, char **argv) {
	int w = -1, h = -1, x0 = w / 2, y0 = h / 2 + 10, ttl = 400;
	double diff = 1, rate = 1, windagl = -M_PI / 2, windvel = 1, a = 1, dt = 1, resolution = 20.;
	string mode = "simulate";
	bool fast = false, draw = true, profiling = false;
	long 	*looptime	= new long[MAX_ITERATIONS],
			*updatetime = new long[MAX_ITERATIONS],
			*optimtime 	= new long[MAX_ITERATIONS],
			*drawtime 	= new long[MAX_ITERATIONS];

	po::options_description desc("All options");
	desc.add_options()
		("help,h", "Displays help message")
		("mode,m", po::value<string>(&mode), "Sets the mode (simulate or concentration)")
		("size,s", po::value<string>(), "Sets the grid dimentions (format w:h, in squares), required")
		("source,S", po::value<string>()->default_value("random"), "Sets the gas source position (format x:y or \"random\")")
		("deltat,t", po::value<double>(&dt)->default_value(dt), "Sets the delta-t of the search, i.e how much time is spent on each step. (in s)")
		("ttl,T", po::value<int>(&ttl)->default_value(ttl), "Sets the particles's expected lifetime (in s)")
		("diffusivity,D", po::value<double>(&diff)->default_value(diff), "Sets the diffusivity of the source (in m²/s)")
		("rate,R", po::value<double>(&rate)->default_value(rate), "Sets the emission rate (in /s, i.e Hz")
		("windangle,A", po::value<double>(&windagl)->default_value(windagl*180/M_PI), "Sets the wind blowing direction (in deg, 0 is EAST, 90 is NORTH)")
		("windvelocity,V", po::value<double>(&windvel)->default_value(windvel), "Sets the velicity of the wind (in m/s)")
		("sensorradius,a", po::value<double>(&a)->default_value(a), "Sets the sensor radius (in m)")
		("resolution,r", po::value<double>(&resolution)->default_value(resolution), "Sets how much grid squares one m equals to (in dots per meter)")
		("fast,f", "Sets the simulator to not use true Infotaxis (Minimize entropy) but instead go where detection is the most probable")
		("nodraw,nd", "Disable frame-per-frame image output. Mean stationary field will still be written.")
#ifndef INFOTAXIS_NO_PROFILING
		("profiling,p", "Enable profiling. Will output a summary of system time taken by respective functions");
		double timeratio;
		{
			struct timespec spec;
			clock_getres(CLOCK_PROCESS_CPUTIME_ID, &spec);
			timeratio = spec.tv_sec + (spec.tv_nsec / 1000000000.);
		}
#else
		;
#endif
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << desc;
		return 0;
	}

	windagl *= M_PI/180; // Back to radians

	if (vm.count("size")) {
		char dummy;
		stringstream(vm["size"].as<string>()) >> w >> dummy >> h;
		if (w < 1 || h < 1) {
			cerr << "Error : Dimensions must be strictly positive" << endl;
			return -1;
		}
	} else {
		cerr << "Argument --size is required" << endl;
		return -1;
	}

	if (vm["source"].as<string>() == "random") {
		srand(time(NULL));
		x0 = rand() % w;
		y0 = rand() % h;
	} else {
		char dummy;
		stringstream(vm["source"].as<string>()) >> x0 >> dummy >> y0;
	}
	if (x0 < 0 || y0 < 0 || x0 >= w || y0 >= h) {
		cerr << "Error : Source coordinates out of bounds" << endl;
		return -1;
	}
	draw = vm.count("nodraw") == 0;
	fast = vm.count("fast") > 0;
#ifdef INFOTAXIS_NO_PROFILING
	if (vm.count("profiling") > 0) {
		cerr << "Enable is disabled. time.h was not found or INFOTAXIS_NO_PROFILING was passed during compilation" << endl;
		return -1;
	}
#else
	profiling = vm.count("profiling") > 0;
#endif

	InfotaxisGrid grid(w, h, diff, rate, windvel, windagl, ttl, a, resolution, !fast);
	struct simul simul = {w / 4, h / 4, 0, &grid};

	if (mode == "simulate") {
		char *filename = new char[64];
		int cnt = 0;
		std::default_random_engine gen;
		gen.seed(time(NULL));

		do {
			PROFILE(looptime[cnt] = (long) -clock();)
			double mean = grid.encounterRate(simul.curx, simul.cury, x0, y0) * dt;
			cout << cnt << " : " << simul.curx << ":" << simul.cury << " => mean " << mean;
			poisson_distribution<int> pdist(mean);
			gen(); // Creating a variate_generator with gen seeds it from gen's value without modifying gen.
			//Thus we need to iterate gen to avoid reseeding the generator with the same value every iteration
			int detects = pdist(gen);

			PROFILE(updatetime[cnt] = (long)-clock();)
			grid.updateProbas(simul.curx, simul.cury, detects, simul.time);
			PROFILE(updatetime[cnt] += clock();)

			cout << ", detections : " << detects << ", entropy : " << grid.entropy() << endl;
			simul.time += dt;

			PROFILE(optimtime[cnt] = (long) -clock();)
			Direction optimal = grid.getOptimalMove(simul.curx, simul.cury, dt);
			PROFILE(optimtime[cnt] += clock();)

			//cout << ", optimal move : "<< DIRECTIONS[optimal+1] << "("<< optimal <<")" << endl;
			sprintf(filename, "pictures/iteration%04d.png", cnt);
			if (draw) {
				PROFILE(drawtime[cnt] = (long)-clock());
				writeProbas(filename, simul, 5, NULL);
				PROFILE(drawtime[cnt] += clock();)
			}
			simul.backtrace.push_front({simul.curx, simul.cury, detects});
			go_to(optimal, simul.curx, simul.cury);
			PROFILE(looptime[cnt] += clock();)
			PROFILE(printf("Time : loop=%5.5fms, update=%5.5fms, opmital find=%5.5fms, drawtime=%-10sms\n", MSECONDS(looptime[cnt]), MSECONDS(updatetime[cnt]), MSECONDS(optimtime[cnt]), (draw ? to_string(MSECONDS(drawtime[cnt])) : "N/A").c_str());)
		} while ((simul.curx != x0 || simul.cury != y0) && cnt++ < MAX_ITERATIONS);
		sprintf(filename, "pictures/iteration%04d.png", cnt);
		if (draw)
			writeProbas(filename, simul, 5, NULL, NULL);

		if (simul.curx == x0 && simul.cury == y0)
			cout << "Source found after " << cnt << " iterations at " << x0 << ":" << y0 << " !" << endl;
		else
			cout << "Source not found after " << cnt << " iterations" << endl;


		if (profiling) {
			for (int i = 1; i < cnt; i ++) {
				looptime[0] += looptime[i];
				updatetime[0] += updatetime[i];
				optimtime[0] += optimtime[i];
				drawtime[0] += drawtime[i];
			}
			printf("Mean time : \n\tloop = % 5.5fms\n\tupdate = % 5.5fms\n\topmital find = % 5.5fms\n\tdrawtime = %-10sms\n", MSECONDS(looptime[0]/cnt), MSECONDS(updatetime[0]/cnt), MSECONDS(optimtime[0]/cnt), (draw ? to_string(MSECONDS(drawtime[0]/cnt)) : "N/A").c_str());
		}
	} else if (vm["mode"].as<string>() != "concentration") {
		cerr << "Unrecognized --mode : "<< vm["mode"].as<string>() << endl;
		return -1;
	}

	writeProbas("pictures/meanfield.png", simul, 5, &x0, &y0);



	cout << "Options : --mode="<< mode << " --size="<< w <<":"<< h <<" --source="<< x0 <<":"<< y0 <<" --deltat="<< dt 
			<<" --ttl="<< ttl <<" --diffusivity="<< diff <<" --rate="<< rate <<" --windangle="<< windagl 
			<<" --windvelocity="<< windvel <<" --sensorradius="<< a <<" --resolution="<< resolution << (fast ? " --fast" : "") << (draw ? "" : " --nodraw") <<endl;

}

void writeProbas(string filename, struct simul &simul, const int ratio, const int *x0, const int *y0) {
	const int imgw = simul.grid->getWidth() * ratio, imgh = simul.grid->getHeight() * ratio;
	image<rgb_pixel> image((size_t) imgw, (size_t) imgh);

	cout << filename << " ";

	if (x0 == NULL || y0 == NULL)
		simul.grid->writeProbabilityField(image, ratio);
	else
		simul.grid->writeMeanStationaryField(image, *x0, *y0, ratio);

	if (simul.backtrace.size() > 0) {
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
	}
	image.write(filename);
}