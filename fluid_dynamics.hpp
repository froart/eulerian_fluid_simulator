#ifndef FLUID_DYNAMICS_HPP_
#define FLUID_DYNAMICS_HPP_

#include <vector>

using namespace std;

class Fluid {
	private:
		float cell_size_;
		int nx_;
		int ny_;
		float dt_;
		int it_;
		float dens_;
		float over_relaxation_;
		vector<float> v_;
		vector<float> u_;
		vector<float> v1_;
		vector<float> u1_;
		vector<float> p_;
		vector<float> s_;
		vector<float> m_;
		vector<float> m1_;
//		float* m_;
//		float* m1_;
		void project();
		void extrapolate();
		void advect_velocity();
		void advect_smoke();
		float sample_field(float, float, int);
	public:
		Fluid(float*, float, int, int, float, float, int);
		void addSmoke(int, int, float);
		void addWind(int, int, float, float);
		float* evaluate();
};

#endif // FLUID_DYNAMICS_HPP_
