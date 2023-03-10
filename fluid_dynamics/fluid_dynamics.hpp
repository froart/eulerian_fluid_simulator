#pragma once

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
		void project_fluid();
		void extrapolate();
		void advect_velocity();
		void advect_smoke();
		float sample_field(float, float, int, vector<float>&);
	public:
		Fluid(float, int, int, float, float, int);
		float* setImage();
		void addSmoke(int, int, float);
		void addWind(int, int, float, float);
		void evaluate();
		void clearImage();
};

