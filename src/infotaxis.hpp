#ifndef INFOTAXIS_HPP_INCLUDED
#define INFOTAXIS_HPP_INCLUDED

#include <future>
#include <list>
#ifndef INFOTAXIS_NOPNG
#include <png++/png.hpp>
#endif

namespace infotaxis {

	typedef struct direction {
		double dx, dy;
		void go_to(double &x, double &y) const {
			x += this->dx;
			y += this->dy;
		}
	} Direction;

	typedef struct simultrace {
		double x, y;
		int detections;
		double dt;
	} Trace;

	enum class ResizeMethod {
		FULL_RECALCULATE, MIXED, RESCALE, FLAT_NEW
	};

	class InfotaxisGrid
	{
	public:
		InfotaxisGrid(int width, int height, double diff, double rate, double windvel, double windang, double part_lifetime, double sensor_radius, double resolution_, bool trueInfotaxis = true);
		InfotaxisGrid(const InfotaxisGrid &grid); // Clone constructor
		~InfotaxisGrid();
		double *operator[](int x);
		double concentration(double x, double y, double x0, double y0);
		double encounterRate(double x, double y, double x0, double y0);
		double expectedEncounterRate(double x, double y);
		double entropy();
		void updateProbas(double x, double y, int n, double t);
		Direction getOptimalMove(double x, double y, double dt, const std::vector<Direction> &directions);
		double getMoveValue(double x, double y, double dt, const Direction &direction);
		double deltaEntropy(double i, double j, double dt);
		int getWidth() const { return width_; };
		int getHeight() const { return height_; };
		int getResolution() const {return resolution_; };
		InfotaxisGrid *resize(int width, int height, int xoff, int yoff, ResizeMethod method = ResizeMethod::MIXED, int rescale_dim = 32);
		void toCSV(std::ofstream& file, std::string separator = ",");
		void toCSV(std::string filename, std::string separator = ","){
			auto file = std::ofstream(filename, std::ios_base::out);
			toCSV(file, separator);
			file.close();}
#ifdef PNGPP_PNG_HPP_INCLUDED
		void writeProbabilityField(png::image<png::rgb_pixel> &image, int ratio);
		void writeMeanStationaryField(png::image<png::rgb_pixel> &image, double x0, double y0, int ratio);
#endif
		static double poisson(double mean, int k);
		const std::vector<Direction> getDefaultDirs() const;
	private:
		InfotaxisGrid(const InfotaxisGrid &grid, int width, int height);
		double *grid_,
				diff_, rate_, windvel_, windang_,
				part_lifetime_, alpha_, lambda_,
				sensor_radius_, resolution_,
				entropyCache;
		int width_, height_;
		bool trueInfotaxis_;
		std::vector<Trace> backtrace_;
		std::vector<Direction> DEFAULTDIRS_;
		void drawArray(png::image<png::rgb_pixel> &image, const double *array, const int w, const int h, const int ratio, const bool absolute);
	};
}

#endif
