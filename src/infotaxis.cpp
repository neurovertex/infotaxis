#include <string.h>
#include <cmath>
#include <iostream>
#include <functional>
#include <thread>
#include <vector>
#include <assert.h>

#include "infotaxis.hpp"

using namespace std;
using namespace png;

namespace infotaxis {

	InfotaxisGrid::InfotaxisGrid(int width, int height, double diff, double rate, double windvel, double windang, double part_lifetime, double sensor_radius, double resolution, bool trueInfotaxis)
	{
		this->width_ = width;
		this->height_ = height;
		this->grid_ = new double[height*width];
		this->diff_ = diff;
		this->rate_ = rate;
		this->windvel_ = windvel;
		this->windang_ = windang;
		this->part_lifetime_ = part_lifetime;
		this->sensor_radius_ = sensor_radius;
		this->resolution_ = resolution; // Grid resolution : cells per meter
		this->trueInfotaxis_ = trueInfotaxis;
		this->entropyCache = -1;
		this->DEFAULTDIRS_.push_back({ 1, 0});
		this->DEFAULTDIRS_.push_back({ 0, 1});
		this->DEFAULTDIRS_.push_back({-1, 0});
		this->DEFAULTDIRS_.push_back({ 0,-1});

		this->alpha_ = rate/(2*M_PI*diff);
		this->lambda_ = sqrt(diff*part_lifetime / (1 + windvel*windvel*part_lifetime/(4*diff)));

		assert(lambda_ > sensor_radius_); // Encounter rate equation won't work if ln(lambda_/sensor_radius_) <= 0

		double val = 1.0/(height*width);
		for (int i = 0; i < height*width; i ++) {
			grid_[i] = val;
		}
	}

	const std::vector<Direction> InfotaxisGrid::getDefaultDirs() const {
		return this->DEFAULTDIRS_;
	}

	InfotaxisGrid::InfotaxisGrid(const InfotaxisGrid &grid)
	{
		this->width_ = grid.width_;
		this->height_ = grid.height_;
		this->grid_ = new double[height_*width_];
		this->diff_ = grid.diff_;
		this->rate_ = grid.rate_;
		this->windvel_ = grid.windvel_;
		this->windang_ = grid.windang_;
		this->part_lifetime_ = grid.part_lifetime_;
		this->sensor_radius_ = grid.sensor_radius_;
		this->resolution_ = grid.resolution_;
		this->alpha_ = grid.alpha_;
		this->lambda_ = grid.lambda_;

		memcpy(this->grid_, grid.grid_, sizeof(double) * width_ * height_);
	}

	InfotaxisGrid::~InfotaxisGrid() {
		delete[] grid_;
	}

	double *InfotaxisGrid::operator[](int y)
	{
		return grid_+(width_ * y);
	}

	double InfotaxisGrid::concentration(double x, double y, double x0, double y0)
	{
		double dist = max(0.001, hypot(x-x0, y-y0)) / resolution_;
		double val = rate_/(4*M_PI*diff_*dist) * exp((sin(windang_)*(y-y0) + cos(windang_)*(x-x0))/resolution_*windvel_/(2*diff_) - dist/lambda_);
		return val;
	}

	double InfotaxisGrid::encounterRate(double x, double y, double x0, double y0)
	{
		double c = concentration(x, y, x0, y0);
		return 2*M_PI*diff_*c / log(lambda_/sensor_radius_);
	}

	double InfotaxisGrid::expectedEncounterRate(double x, double y)
	{
		double sum = 0;
		for (int i = 0; i < width_*height_; i ++) {
			double er = encounterRate(x, y, i % width_, i / width_);
			if (std::isfinite(er))
				sum += grid_[i] * er;
			else
				cerr << "NaN encounterRate at "<< x <<":"<< y <<"|"<< i%width_ <<":"<< i/width_<< endl;
		}
		return sum;
	}

	void InfotaxisGrid::updateProbas(double x, double y, int n, double dt)
	{
		if (!isfinite(x) || !isfinite(y) || !isfinite(dt) || n<0) {
			printf("updateProbas(x=%f, y=%f, n=%d, dt=%f)\n", x, y, n, dt);
			cerr << "Error : invalid parameters" << endl;
			throw "Invalid parameters";
		}
		double sum = 0;
		entropyCache = -1; // Invalidate cache
		for (int j = 0; j < height_; j ++)
			for (int i = 0; i < width_; i ++)
				sum += grid_[i + width_*j] *= pow(encounterRate(x, y, i, j)*dt, n) * exp(- encounterRate(x, y, i, j)*dt);


		sum = 1/sum;

		for (int j = 0; j < height_; j ++) for (int i = 0; i < width_; i ++)
			grid_[i + width_*j] *= sum;

		backtrace_.push_back({x, y, n});
	}

	double InfotaxisGrid::entropy()
	{
		if (entropyCache < 0) {
			double sum = 0;

			for (int i = 0; i < height_*width_; i ++)
				sum -= (grid_[i] > 0) ? grid_[i] * log(grid_[i]) : 0; // 0*log(0) => NaN, but lim[x->0](x*log(x)) = 0 and that's what matters
			entropyCache = sum;
		}
		return entropyCache;
	}

	double InfotaxisGrid::poisson(const double mean, const int k)
	{
		return (k == 0) ? exp(-mean) : poisson(mean, k-1) * mean / k;
	}

	double InfotaxisGrid::deltaEntropy(const double i, const double j, const double dt)
	{
		InfotaxisGrid *grid = this;
		double prob = grid->grid_[grid->width_ * (int)j + (int)i], delta = 0, cumul = 0, p,
				ear = grid->expectedEncounterRate(i, j), mean = ear * dt, entropy = grid->entropy();
		InfotaxisGrid newGrid = InfotaxisGrid(*grid);
		for (int k = max((int)mean - 10, 0); cumul < 0.9999 && k < mean + 10; k ++) {
			cumul += p = InfotaxisGrid::poisson(mean, k);
			memcpy(newGrid.grid_, grid->grid_, sizeof(double) * grid->height_ * grid->width_);
			newGrid.updateProbas(i, j, k, dt);
			if (!isfinite(p)) {
				cerr << "Non-finite poisson : "<< i <<":"<< j <<" : "<< p <<"("<< mean <<", "<< k <<")" << endl;
				assert(false);
			}
			double e = newGrid.entropy();
			if (!isfinite(e)) {
				cerr << "Non-finite entropy : "<< i <<":"<< j <<" : "<< e << endl;
				assert(false);
			}
			delta += p * (e - entropy);
		}
		delta =  prob * entropy - (1-prob)* delta;
		return delta;
	}

	Direction InfotaxisGrid::getOptimalMove(const double x, const double y, const double dt, const vector<Direction> &directions)
	{
		double best = -1;
		Direction bestDir = directions[0];
		vector<pair<Direction, future<double>>> promises;
		for (Direction d : directions) {
			double i = x, j = y;
			d.go_to(i, j);
			if (i >= 0 && i < width_ && j >= 0 && j < height_) {
				if (trueInfotaxis_)
					promises.push_back(pair<Direction, future<double>>(d,
						async(launch::async, &InfotaxisGrid::deltaEntropy, this,  i, j, dt)));
				else
					promises.push_back(pair<Direction, future<double>>(d,
						async(launch::async, &InfotaxisGrid::expectedEncounterRate, this,  i, j)));
			}
		}

		for (auto&& p : promises) {
			double val = get<1>(p).get();
			assert(isfinite(val));
			if (val > best) {
				best = val;
				bestDir = get<0>(p);
			}
		}

		if (best == 1)
			cerr << "Error : no direction results in information gain (not supposed to happen)" << endl; // Yeah sure it's not, well I'm gonna check anyway
		return bestDir;
	}

	double InfotaxisGrid::getMoveValue(double x, double y, double dt, const Direction &direction) {
		direction.go_to(x, y);
		if (x >= 0 && y >= 0 && x < width_ && y < height_)
			return deltaEntropy(x, y, dt);
		else
			return 0;
	}

#ifdef PNGPP_PNG_HPP_INCLUDED
    using namespace png;

	rgb_pixel getColour(double absv, double relv)
	{
		double value = min((int)(relv*32)*8, 255);
		byte b = (byte)(value * cos(absv*M_PI/2)),
				r = (byte)(value * sin(absv*M_PI/2)),
				g = b;
		return rgb_pixel(r, g, b);
	}

	void InfotaxisGrid::drawArray(image<rgb_pixel> &image, const double *array, const int w, const int h, const int ratio, const bool absolute = true)
	{
		const int imgw = w * ratio, imgh = h * ratio;
		double mi = INFINITY, ma = -INFINITY;
		double *values = new double[w*h];
		for (int i = 0; i < h; i ++)
			for (int j = 0; j < w; j ++) {
				double raw = array[i*w + j], val = (raw);
				if (std::isinf(raw))
					cerr << "Infinite : "<< j << ":"<< i << endl;
				else if (std::isnan(raw))
					cerr << "NaN : "<< j << ":"<< i << endl;
				if (!std::isnan(val) && !std::isinf(val)) {
					mi = min(mi, val);
					ma = max(ma, val);
				} else
					val = min(max(val, mi), ma);
				values[i*w + j] = val;
			}

		double scale = (ma != mi) ? ma-mi : 1;

		for (int i = 0; i < imgh; i ++)
			for (int j = 0; j < imgw; j ++) {
				image[imgh-1-i][j] = getColour(
						absolute ? values[i/ratio*w + j/ratio] : 0,
						(values[i/ratio*w + j/ratio]-mi)/scale);
			}
		delete[] values;

		if (backtrace_.size() > 0) {
			for (Trace &t : backtrace_) {
				rgb_pixel color = (t.detections == 0) ? rgb_pixel(64, 128, 128) : rgb_pixel(
						(byte) (255 - (128 / t.detections)), 64, 64);
				for (int j = 1; j < max(ratio - 1, 2); j++)
					for (int i = 1; i < max(ratio - 1, 2); i++) {
						image[imgh - 1 - ((int)(t.y) * ratio + (int)j)][t.x * ratio + i] = color;
					}
			}
			Trace last = backtrace_.back();
			int x = (int)(last.x + .5) * ratio;
			int y = (int)(last.y + .5) * ratio;
			for (int i = 0; i <= 10; i++) {
				if (imgh - y + 4 - i >= 0 && imgh - y + 4 - i < imgh) {
					if (x - 5 + i >= 0 && x - 5 + i < imgw)
						image[imgh - y + 4 - i][x - 5 + i] = rgb_pixel(255, 128, 0);

					if (x + 5 - i >= 0 && x + 5 - i < imgw)
						image[imgh - y + 4 - i][x + 5 - i] = rgb_pixel(255, 128, 0);
				}
			}
		}
	}

	void InfotaxisGrid::writeProbabilityField(png::image<png::rgb_pixel> &image, const int ratio)
	{
		this->drawArray(image, grid_, width_, height_, ratio);
	}
	void InfotaxisGrid::writeMeanStationaryField(png::image<png::rgb_pixel> &image, const double x0, const double y0, const int ratio)
	{
		double *array = new double[width_*height_];
		for (int y = 0; y < height_; y ++)
			for (int x = 0; x < width_; x ++)
				array[y * width_ + x] = log(concentration(x, y, x0, y0));

		this->drawArray(image, array, width_, height_, ratio, false);
	}

#endif
}
