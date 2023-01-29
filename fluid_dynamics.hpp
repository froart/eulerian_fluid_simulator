#ifndef FLUID_DYNAMICS_HPP_
#define FLUID_DYNAMICS_HPP_

#include <memory>

using namespace std;

class Fluid {
	private:
		int w_; 
		int h_;	
		float dt_;
		int it_;
		float* d_;
		unique_ptr<float[]> v_;
		unique_ptr<float[]> u_;
		unique_ptr<float[]> v1_;
		unique_ptr<float[]> u1_;
		unique_ptr<float[]> p_;
		unique_ptr<float[]> s_;
		unique_ptr<float[]> m_;
	public:
		Fluid(float*, int, int, float, int);
		~Fluid();
		void addSmoke(int, int, float);
		void addWind(int, int, float, float);
		void evaluate();
};

#endif // FLUID_DYNAMICS_HPP_
