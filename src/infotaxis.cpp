#include <string.h>
#include <cmath>
#include <iostream>
#include <assert.h>
#include <boost/math/distributions/poisson.hpp>
//#include <boost/math/special_functions/bessel.hpp>

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
		for (int j = 0; j < height_; j ++) for (int i = 0; i < width_; i ++)
			sum += grid_[i + width_*j] *= pow(encounterRate(x, y, i, j)*dt, n) * exp(- encounterRate(x, y, i, j)*dt);


		for (int j = 0; j < height_; j ++) for (int i = 0; i < width_; i ++)
			grid_[i + width_*j] /= sum;

		setLastTime(t);
	}

	double InfotaxisGrid::entropy()
	{
		double sum = 0;

		for (int i = 0; i < height_*width_; i ++)
			sum -= (grid_[i] > 0) ? grid_[i] * log(grid_[i]) : 0; // 0*log(0) => NaN, but lim[x->0](x*log(x)) = 0 and that's what matters

		return sum;
	}

	Direction InfotaxisGrid::getOptimalMove(const int x, const int y, const double dt)
	{
		const double entropy = this->entropy();
		double best = 1;
		Direction bestDir = STAY;

		for (int d = STAY; d <= SOUTH; d ++) {
			int i = x, j = y;
			go_to((Direction)d, i, j);
			if (i >= 0 && i < width_ && j >= 0 && j < height_) {
				double prob = grid_[width_ * j + i], delta = 0, cumul = 0, p = 1, ear = expectedEncounterRate(i, j);
				boost::math::poisson_distribution<double> distr(ear * dt);
				InfotaxisGrid newGrid = InfotaxisGrid(*this);
				//cout << DIRECTIONS[d+1] <<" \t: "<< distr.mean() << endl;
				for (int k = max(0, (int)distr.mean()-5); k < max(0, (int)distr.mean()-5)+10 && p > 0.01; k ++) {
					memcpy(newGrid.grid_, this->grid_, sizeof(double) * height_ * width_);
					newGrid.setLastTime(last_time_);
					newGrid.updateProbas(i, j, k, last_time_ + dt);
					p = pdf(distr, k);
					double e = newGrid.entropy();
					delta += p * (e - entropy);
					//cout << "\tdelta : "<< delta << "\t(p=" << p << "\t,e="<< e <<")" << endl;
				}
				delta = (1-prob)* delta - prob * entropy;
				//cout << "\tFinal delta : " << delta << endl;
				if (delta < best) {
					best = delta;
					bestDir = (Direction) d;
				}
			}
		}
		cout << "Best dS : "<< DIRECTIONS[bestDir +1] << " : "<< best << endl;
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
}