#include "infotaxis.hpp"
#include <stdexcept>
#include <algorithm>
#include <list>
#include <boost/program_options.hpp>

#ifdef INFOTAXIS_NO_PROFILING
#define PROFILE(cmd)
#define MSECONDS(time) (double)(.0)
#else
#include <time.h>
#define PROFILE(cmd) if (profiling) cmd;
#define MSECONDS(time) (((time) * timeratio * 1000.))
#endif

#define MAX_ITERATIONS 2048

using namespace std;
using namespace png;
using namespace infotaxis;
namespace po = boost::program_options;

struct simul {
	double curx, cury;
	InfotaxisGrid *grid;
	list<Trace> backtrace;
};

void writeProbas(string filename, InfotaxisGrid &grid, const int ratio, const int *x0 = nullptr, const int *y0 = nullptr);

void simulate(struct simul &simul, int x0, int y0, double dt, bool profiling, bool draw, ofstream *logfile, bool quiet);

static const int RATIO = 2;

int main(int argc, char **argv) {
	int w = -1, h = -1, x0, y0, ttl = 2500;
	double diff = 1, rate = 1, windagl = -M_PI / 2, windvel = 1, a = .1, dt = 1, resolution = 20., startx, starty;
	ofstream *logfile = NULL;
	string meanfield;
	bool nosimulate = false, fast = false, draw = true, profiling = false, quiet;

	po::options_description desc("All options");
	desc.add_options()
		("help,h", "Displays help message")
		("nosimulate", "Disables simulation")
		("meanfield,m", po::value<string>(&meanfield)->default_value(""), "Prints the mean stationary field of the gas concentration to <file>")
		("size", po::value<string>(), "Sets the grid dimentions (format w:h, in squares), required")
		("start,s", po::value<string>()->default_value("random"), "Sets the robot initial position (format x:y or \"random\")")
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
		("nodraw,n", "Disable frame-per-frame image output. Mean stationary field will still be written.")
		("quiet,q", "Disable all normal output, only print final iteration count after simulation.")
		("log,l", po::value<string>(), "Logs position and entropy")
#ifndef INFOTAXIS_NO_PROFILING
		("profiling,p", "Enable profiling. Will output a summary of system time taken by respective functions");
#else
		;
#endif
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);

	if (vm.count("help") > 0) {
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

	if (vm["start"].as<string>() == "random") {
		srand(time(NULL));
		startx = rand() % w;
		starty = rand() % h;
	} else {
		char dummy;
		stringstream(vm["start"].as<string>()) >> startx >> dummy >> starty;
	}
	if (startx < 0 || starty < 0 || startx >= w || starty >= h) {
		cerr << "Error : Initial coordinates out of bounds" << endl;
		return -1;
	}

	nosimulate = vm.count("nosimulate") > 0;
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
	logfile = (vm.count("log") > 0) ? new ofstream(vm["log"].as<string>()) : NULL;
	quiet = vm.count("quiet") > 0;

	if (nosimulate && meanfield == "") {
		cerr << "Simulation and meanfield output disabled. Nothing to do" << endl;
		exit(-1);
	}

	InfotaxisGrid grid(w, h, diff, rate, windvel, windagl, ttl, a, resolution, !fast);
	struct simul simul = {startx, starty, &grid};

	if (logfile != NULL)
		*logfile << w <<" "<< h <<" "<< diff <<" "<< rate <<" "<< windvel <<" "<< windagl <<" "<< ttl
			<<" "<< a <<" "<< resolution <<" "<< fast << dt << endl << endl;

	if (!nosimulate)
		try {
			simulate(simul, x0, y0, dt, profiling, draw, logfile, quiet);
		} catch (const exception &e) {
			cerr << "Exception : "<< e.what() << endl;
			return 1;
		}

	if (logfile != NULL) {
		*logfile << flush;
		logfile->close();
	}

	if (draw && !nosimulate && meanfield == "")
		writeProbas("pictures/meanfield.png", *simul.grid, RATIO, &x0, &y0);
	else if (meanfield != "")
		writeProbas(meanfield, *simul.grid, RATIO, &x0, &y0);

}

void simulate(struct simul &simul, int x0, int y0, double dt, bool profiling, bool draw, ofstream *logfile, bool quiet) {
	char *filename = new char[64];
	long 	*looptime = new long[4 * MAX_ITERATIONS], // Slicing the buffer in 4
			*updatetime = looptime + MAX_ITERATIONS,
			*optimtime = updatetime + MAX_ITERATIONS,
			*drawtime = optimtime + MAX_ITERATIONS;
	InfotaxisGrid &grid = *simul.grid;
	int cnt = 0;
	std::default_random_engine gen;
	gen.seed(time(NULL));
#ifndef INFOTAXIS_NO_PROFILING
	double timeratio;
	{
		struct timespec spec;
		clock_getres(CLOCK_PROCESS_CPUTIME_ID, &spec);
		timeratio = spec.tv_sec + (spec.tv_nsec / 1000000000.);
	}
#endif

	InfotaxisGrid *resized = nullptr;
	do {
		PROFILE(looptime[cnt] = (long) -clock();)
		double mean = grid.encounterRate(simul.curx, simul.cury, x0, y0) * dt, entropy = grid.entropy();
		poisson_distribution<int> pdist(mean);
		gen(); // Creating a variate_generator with gen seeds it from gen's value without modifying gen.
		//Thus we need to iterate gen to avoid reseeding the generator with the same value every iteration
		int detects = pdist(gen);

		if (!quiet)
			cout << cnt << " : " << simul.curx << ":" << simul.cury << " => mean " << mean <<
				", detections : " << detects << ", entropy : " << entropy << endl;

		if (logfile != NULL)
			*logfile << simul.curx <<" "<< simul.cury <<" "<< entropy <<" "<< detects << endl;

		PROFILE(updatetime[cnt] = (long)-clock();)
		grid.updateProbas(simul.curx, simul.cury, detects, dt);
		PROFILE(updatetime[cnt] += clock();)

		PROFILE(optimtime[cnt] = (long) -clock();)
		Direction optimal = grid.getOptimalMove(simul.curx, simul.cury, dt, grid.getDefaultDirs());
		PROFILE(optimtime[cnt] += clock();)

		sprintf(filename, "pictures/iteration%04d.png", cnt);
		if (draw) {
			PROFILE(drawtime[cnt] = (long)-clock());
			writeProbas(filename, *simul.grid, RATIO);
			PROFILE(drawtime[cnt] += clock();)
		}
		simul.backtrace.push_front({simul.curx, simul.cury, detects});
		optimal.go_to(simul.curx, simul.cury);
		PROFILE(looptime[cnt] += clock();)
		PROFILE(printf("Time : loop=%5.5fms, update=%5.5fms, opmital find=%5.5fms, drawtime=%-10sms\n", MSECONDS(looptime[cnt]), MSECONDS(updatetime[cnt]), MSECONDS(optimtime[cnt]), (draw ? to_string(MSECONDS(drawtime[cnt])) : "N/A").c_str());)
	} while (cnt++ < MAX_ITERATIONS && (simul.curx != x0 || simul.cury != y0));
	sprintf(filename, "pictures/iteration%04d.png", cnt);
	if (draw)
		writeProbas(filename, *simul.grid, RATIO);

	if (resized != nullptr)
		delete resized;

	if (!quiet) {
		if (simul.curx == x0 && simul.cury == y0)
			cout << "Source found after " << cnt << " iterations at " << x0 << ":" << y0 << " !" << endl;
		else
			cout << "Source not found after " << cnt << " iterations" << endl;
	} else {
		cout << cnt << flush;
	}


	if (profiling) {
		for (int i = 1; i < cnt; i ++) {
			looptime[0] += looptime[i];
			updatetime[0] += updatetime[i];
			optimtime[0] += optimtime[i];
			drawtime[0] += drawtime[i];
		}
		printf("Mean time : \n\tloop = % 5.5fms\n\tupdate = % 5.5fms\n\topmital find = % 5.5fms\n\tdrawtime = %-10sms\n", MSECONDS(looptime[0]/cnt), MSECONDS(updatetime[0]/cnt), MSECONDS(optimtime[0]/cnt), (draw ? to_string(MSECONDS(drawtime[0]/cnt)) : "N/A").c_str());
	}
}

void writeProbas(string filename, InfotaxisGrid &grid, const int ratio, const int *x0, const int *y0) {
	const int imgw = grid.getWidth() * ratio, imgh = grid.getHeight() * ratio;
	image<rgb_pixel> image((size_t) imgw, (size_t) imgh);
	if (x0 == NULL || y0 == NULL)
		grid.writeProbabilityField(image, ratio);
	else
		grid.writeMeanStationaryField(image, *x0, *y0, ratio);

	image.write(filename);
}
