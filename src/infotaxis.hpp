#ifndef INFOTAXIS_HPP_INCLUDED
#define INFOTAXIS_HPP_INCLUDED

#include <future>
#ifndef INFOTAXIS_NOPNG
#include <png++/png.hpp>
#endif

namespace infotaxis {

	static const char * const DIRECTIONS[] = {"STAY", "EAST", "NORTH", "WEST", "SOUTH"};

	typedef enum direction {
		STAY=-1,
		EAST=0,
		NORTH=1,
		WEST=2,
		SOUTH=3
	} Direction;

	void go_to(Direction d, int &x, int &y);

	class InfotaxisGrid
	{
	public:
		InfotaxisGrid(int width, int height, double diff, double rate, double windvel, double windang, double part_lifetime, double sensor_radius, double resolution_, bool trueInfotaxis = true);
		InfotaxisGrid(const InfotaxisGrid &grid); // Clone constructor
		~InfotaxisGrid();
		double *operator[](int x);
		double concentration(int x, int y, int x0, int y0);
		double encounterRate(int x, int y, int x0, int y0);
		double expectedEncounterRate(int x, int y);
		double entropy();
		void updateProbas(int x, int y, int n, double t);
		void setLastTime(double t);
		Direction getOptimalMove(int x, int y, double dt);
		int getWidth() { return width_; };
		int getHeight() { return height_; };
#ifdef PNGPP_PNG_HPP_INCLUDED
		void writeProbabilityField(png::image<png::rgb_pixel> &image, int ratio);
		void writeMeanStationaryField(png::image<png::rgb_pixel> &image, int x0, int y0, int ratio);
#endif
		static double poisson(double mean, int k);
	private:
		double *grid_,
				diff_, rate_, windvel_, windang_, 
				part_lifetime_, alpha_, lambda_, 
				sensor_radius_, last_time_, resolution_,
				entropyCache;
		int width_, height_;
		bool trueInfotaxis_;
		double deltaEntropy(const int i, const int j, const double dt);
	};
}

#endif