namespace infotaxis {

	static const char *DIRECTIONS[] = {"STAY", "EAST", "NORTH", "WEST", "SOUTH"};

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
		InfotaxisGrid(int width, int height, double diff, double rate, double windvel, double windang, double part_lifetime, double sensor_radius, double resolution_);
		InfotaxisGrid(const InfotaxisGrid &grid); // Clone constructor
		double *operator[](int x);
		double concentration(int x, int y, int x0, int y0);
		double encounterRate(int x, int y, int x0, int y0);
		double expectedEncounterRate(int x, int y);
		double entropy();
		void updateProbas(int x, int y, int n, double t);
		void setLastTime(double t);
		Direction getOptimalMove(int x, int y, double dt);
	private:
		double *grid_,
				diff_, rate_, windvel_, windang_, 
				part_lifetime_, alpha_, lambda_, 
				sensor_radius_, last_time_, resolution_;
		int width_, height_;
	};
}