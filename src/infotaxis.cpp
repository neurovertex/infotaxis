#include <string.h>
#include <cmath>
#include <iostream>
#include <functional>
#include <thread>
#include <vector>
#include <assert.h>

#include "infotaxis.hpp"

using namespace std;

namespace infotaxis {

	InfotaxisGrid::InfotaxisGrid(int width, int height, double diff, double rate, double windvel, double windang, double part_lifetime, double sensor_radius, double resolution) 
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
		this->last_time_ = 0;
		this->resolution_ = resolution; // Grid resolution : edge size of its squares, in meters.
		this->entropyCache = -1;

		this->alpha_ = rate/(2*M_PI*diff);
		this->lambda_ = sqrt(diff*part_lifetime / (1 + windvel*windvel*part_lifetime/(4*diff)));

		cout << "alpha : "<< alpha_ << ", lambda : "<< lambda_ << endl;

		assert(lambda_ > sensor_radius_); // Encounter rate equation won't work if ln(lambda_/sensor_radius_) <= 0

		double val = 1.0/(height*width);
		for (int i = 0; i < height*width; i ++) {
			grid_[i] = val;
		}
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
		this->last_time_ = grid.last_time_;
		this->resolution_ = grid.resolution_;
		this->alpha_ = grid.alpha_;
		this->lambda_ = grid.lambda_;

		memcpy(this->grid_, grid.grid_, sizeof(double) * width_ * height_);
	}

	InfotaxisGrid::~InfotaxisGrid() {
		delete grid_;
	}

	void InfotaxisGrid::setLastTime(double t)
	{
		this->last_time_ = t;
	}

	double *InfotaxisGrid::operator[](int y)
	{
		return grid_+(width_ * y);
	}

	double InfotaxisGrid::concentration(int x, int y, int x0, int y0)
	{
		double dist = max(0.001, hypot(x-x0, y-y0)) * resolution_;
		double val = rate_/(4*M_PI*diff_*dist) * exp((sin(windang_)*(y-y0) + cos(windang_)*(x-x0))*resolution_*windvel_/(2*diff_) - dist/lambda_);
		return val;
	}

	double InfotaxisGrid::encounterRate(int x, int y, int x0, int y0)
	{
		double c = concentration(x, y, x0, y0);
		return 2*M_PI*diff_*c / log(lambda_/sensor_radius_);
	}

	double InfotaxisGrid::expectedEncounterRate(int x, int y)
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

	void InfotaxisGrid::updateProbas(int x, int y, int n, double t)
	{
		double sum = 0;
		double dt = t - last_time_;
		entropyCache = -1; // Invalidate cache
		for (int j = 0; j < height_; j ++) for (int i = 0; i < width_; i ++)
			sum += grid_[i + width_*j] *= pow(encounterRate(x, y, i, j)*dt, n) * exp(- encounterRate(x, y, i, j)*dt);


		for (int j = 0; j < height_; j ++) for (int i = 0; i < width_; i ++)
			grid_[i + width_*j] /= sum;

		setLastTime(t);
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

	double InfotaxisGrid::deltaEntropy(InfotaxisGrid *grid, const int i, const int j, const double dt)
	{
		double prob = grid->grid_[grid->width_ * j + i], delta = 0, cumul = 0, p, ear = grid->expectedEncounterRate(i, j),
					mean = ear * dt, entropy = grid->entropy();
		InfotaxisGrid newGrid = InfotaxisGrid(*grid);
		for (int k = max((int)mean - 10, 0); cumul < 0.9999 && k < mean + 10; k ++) {
			cumul += p = InfotaxisGrid::poisson(mean, k);
			memcpy(newGrid.grid_, grid->grid_, sizeof(double) * grid->height_ * grid->width_);
			newGrid.setLastTime(grid->last_time_);
			newGrid.updateProbas(i, j, k, grid->last_time_ + dt);
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
			//cout << "\tdelta : "<< delta << "\t(p=" << p << "\t,e="<< e <<")" << endl;
		}
		delta = (1-prob)* delta - prob * entropy;
		return delta;
	}

	Direction InfotaxisGrid::getOptimalMove(const int x, const int y, const double dt)
	{
		double best = 1;
		Direction bestDir = STAY;
		cout << "dS's : {";
		vector<pair<Direction, future<double>>> promises;
		for (int d = STAY; d <= SOUTH; d ++) {
			int i = x, j = y;
			go_to((Direction)d, i, j);
			if (i >= 0 && i < width_ && j >= 0 && j < height_) {
				promises.push_back(pair<Direction, future<double>>((Direction)d, 
					async(launch::async, deltaEntropy, this,  i, j, dt)));
			}
		}

		for (auto&& p : promises) {
			double val = get<1>(p).get();
			cout << DIRECTIONS[get<0>(p)> +1][0] <<":"<< val << ",";
			if (val < best) {
				best = val;
				bestDir = get<0>(p);
			}
		}

		cout << "Best dS : "<< DIRECTIONS[bestDir +1] << " : "<< best << ")" << endl;
		if (best == 1)
			cerr << "Error : no direction results in information gain (not supposed to happen)" << endl; // Yeah sure it's not, well I'm gonna check anyway
		return bestDir;
	}

	void go_to(Direction d, int &x, int &y)
	{
		int dx, dy;
		if (d == STAY) {
	 		dx = dy = 0;
		} else {
			dx = (1-d%2)*(1-2*(d/2)),
			dy = (d%2)*(1-2*(d/2));
		}
		x += dx;
		y += dy;
	}

#ifdef PNGPP_PNG_HPP_INCLUDED
    using namespace png;

	rgb_pixel getColour(double absv, double relv)
	{
		double value = min((int)(relv*32)*8, 255);
		byte b = (byte)(value * cos(absv*M_PI/2)),
				r = (byte)(value * sin(absv*M_PI/2)),
				g = b;
		//cout << "relv : " << relv << ", absv : " << absv << ", rgb : "<< r <<":"<< g <<":"<< b << endl;
		return rgb_pixel(r, g, b);
	}

	void drawArray(image<rgb_pixel> &image, const double *array, const int w, const int h, const int ratio, const bool absolute = true)
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
				//cout << raw << endl;
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
		//		cout << j << ":" << i <<" =>\t" << values[i*w/ratio + j/ratio] << endl;
				image[imgh-1-i][j] = getColour(
						absolute ? values[i/ratio*w + j/ratio] : 0,
						(values[i/ratio*w + j/ratio]-mi)/scale);
			}
		cout << "min-max: "<< mi <<":"<< ma << endl;
	}

	void InfotaxisGrid::writeProbabilityField(png::image<png::rgb_pixel> &image, const int ratio)
	{
		drawArray(image, grid_, width_, height_, ratio);
	}
	void InfotaxisGrid::writeMeanStationaryField(png::image<png::rgb_pixel> &image, const int x0, const int y0, const int ratio)
	{
		double *array = new double[width_*height_];
		for (int y = 0; y < height_; y ++)
			for (int x = 0; x < width_; x ++)
				array[y * width_ + x] = log(concentration(x, y, x0, y0));

		drawArray(image, array, width_, height_, ratio, false);
	}

#endif
}