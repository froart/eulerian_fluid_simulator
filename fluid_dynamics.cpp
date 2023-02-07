#include "fluid_dynamics.hpp"
#include "opengl_setup.hpp"
#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>

using namespace std;

// DEBUG
#define mark(name) cout << "reached the mark; " << name << endl // DEBUG MARK
#define freeze(flag) if(flag){cout << "simulation frozen! "; for(;;);} 

#define U_FIELD 0
#define V_FIELD 1
#define S_FIELD 2

#define I(x,y) ((x)+(nx_)*(y))

Fluid::Fluid(float cell_size,
						 int nx, 
						 int ny, 
						 float dt,
						 float dens,
						 int iterations) 
						 : cell_size_(cell_size),
							 nx_(nx),
							 ny_(ny),
							 dt_(dt),
							 dens_(dens),
							 it_(iterations) {
	int cell_num = nx_ * ny_;
	vector<float> u(cell_num, 0.0);
	vector<float> v(cell_num, 0.0);
	vector<float> u1(cell_num, 0.0);
	vector<float> v1(cell_num, 0.0);
	vector<float> p(cell_num, 0.0);
	vector<float> s(cell_num, 1.0); // 1.0 for fluid
	vector<float> m(cell_num, 0.0);
	vector<float> m1(cell_num, 0.0);
	for(int j = 0; j < ny_; ++j)
		for(int i = 0; i < nx_; ++i) {
			if(i == 0 || j == 0 || j == ny_-1) // right border is opened for the smoke to flow out
				s[I(i,j)] == 0.0; // 0.0 for obstacle
			}
	u_.swap(u);
	v_.swap(v);
	u1_.swap(u1);
	v1_.swap(v1);
	p_.swap(p);
	s_.swap(s);
	m_.swap(m);
	m1_.swap(m1);
	over_relaxation_ = 1.9;
};

float* Fluid::setImage() {
	return &m_[0];
}
int flag = 0;
void Fluid::addSmoke(int x, int y, float amount) {
	m_[I(x,y)] = amount;
	flag = 1;
}

void Fluid::addWind(int x, int y, float x_amount, float y_amount) {
	u_[I(x,y)] = x_amount;
	v_[I(x,y)] = y_amount;
}

void Fluid::evaluate() {
	fill(p_.begin(), p_.end(), 0.0);
	project();
	extrapolate();
	advect_velocity();
	advect_smoke();
}

void Fluid::clearImage() {
	std::fill(m_.begin(), m_.end() , 0.0);
}

void Fluid::project() { // force imcompressibility
	for(int k = 0; k < it_; k++)
		for(int j = 1; j < ny_-1; j++)
			for(int i = 1; i < nx_-1; i++) {
				if(!s_[I(i,j)]) continue; // if evaluating at an obstacle
				float sx0 = s_[I(i-1,j)];
				float sx1 = s_[I(i+1,j)];
				float sy0 = s_[I(i,j-1)];
				float sy1 = s_[I(i,j+1)];
				float s = sx0 + sx1 + sy0 + sy1;
				if(!s) continue; // if nothing to compute here...
				float div = u_[I(i+1,j)] - u_[I(i,j)] + v_[I(i,j+1)] - v_[I(i,j)];
				// computing the pressure (not necessary for the simulation
				// FIXME  should be minus???
				float p = -(div / s) * over_relaxation_;
				// TODO where does this equation come from???
				p_[I(i,j)] += p * (dens_ * cell_size_) / dt_; 
				// fix the fluid to be imcompressible
				u_[I(i,j)] -= sx0 * p;
				u_[I(i+1,j)] += sx1 * p;
				v_[I(i,j)] -= sy0 * p;
				v_[I(i,j+1)] += sy1 * p;
			}
}

void Fluid::extrapolate() { // enforce border conditions
	for(int i = 0; i < nx_; ++i) {
		u_[I(i,0)] = u_[I(i,1)]; // left wall
		u_[I(i,ny_-1)] = u_[I(i,ny_-2)]; // right wall
	}
	for(int j = 0; j < ny_; ++j) {
		v_[I(0,j)] = v_[I(1,j)]; // bottom wall
		v_[I(nx_-1,j)] = v_[I(nx_-2,j)]; // top wall
	}
}

void Fluid::advect_velocity() {
	u1_ = u_;
	v1_ = v_;
	float h = cell_size_;
	float h2 = cell_size_ * 0.5;
	for(int j = 1; j < ny_; j++)
		for(int i = 1; i < nx_; i++) {
			if(s_[I(i,j)] && s_[I(i,j-1)] && j < ny_-1) { // u-component
				float x = (float) i * h;
				float y = (float) j * h + h2; // u-component is situated at this point
				float u = u_[I(i,j)];	
				float v = (v_[I(i-1,j)] + v_[I(i-1,j+1)] + v_[I(i,j)] + v_[I(i,j+1)]) * 0.25;	// v-component in this case is averaged by 4 values around
				x -= u*dt_;
				y -= v*dt_;
				u = sample_field(x, y, U_FIELD, u_);
				u1_[I(i,j)] = u;
			}
			if(s_[I(i,j)] && s_[I(i-1,j)] && i < nx_-1) { // v-component
				float x = (float) i * h + h2; // v-component is situated at this point
				float y = (float) j * h;
				float u = (u_[I(i,j)] + u_[I(i,j-1)] + u_[I(i+1,j-1)] + u_[I(i+1,j)]) * 0.25;	// u-component in this case is averaged by 4 values around
				float v = v_[I(i,j)];	
				x -= u*dt_;
				y -= v*dt_;
				v = sample_field(x, y, V_FIELD, v_);
				v1_[I(i,j)] = v;
			}
		}
	u_ = u1_;
	v_ = v1_;
}

void Fluid::advect_smoke() {
	m1_= m_;	
	float h = cell_size_;
	float h2 = 0.5 * h;
	for(int j = 1; j < ny_-1; ++j)
		for(int i = 1; i < nx_-1; ++i)
			if(s_[I(i,j)]) {
				float u = (u_[I(i,j)] + u_[I(i+1,j)]) * 0.5;	
				float v = (v_[I(i,j)] + v_[I(i,j+1)]) * 0.5;	
				float x = ((float) i) * h + h2 - dt_ * u;
				float y = ((float) j) * h + h2 - dt_ * v;
				m1_[I(i,j)] = sample_field(x, y, S_FIELD, m_);
			}
	m_= m1_;	
}

float Fluid::sample_field(float x_p, float y_p, int field, vector<float>& vec) {
	float h = cell_size_;
	float h1 = 1.0 / cell_size_;
	float h2 = cell_size_ * 0.5;
	
	float x = fmax(fmin(x_p, nx_*h), h);
	float y = fmax(fmin(y_p, ny_*h), h);

	float dx, dy = 0.0;
	switch(field) {
		case U_FIELD: dy = h2; break;
		case V_FIELD: dx = h2; break;
		case S_FIELD: dx = h2; dy = h2; break;
	}
	// TODO: what is it?
	float x0 = fmin(floor((x-dx)*h1), nx_-1); 	
	float tx = ((x-dx) - x0*h) * h1;
	float x1 = fmin(x0 + 1, nx_-1);

	float y0 = fmin(floor((y-dy)*h1), ny_-1); 	
	float ty = ((y-dy) - y0*h) * h1;
	float y1 = fmin(y0 + 1, ny_-1);

	float sx = 1.0 - tx;
	float sy = 1.0 - ty;

	float val = sx * sy * vec[I(x0,y0)]
						+ tx * sy * vec[I(x1,y0)]
						+ tx * ty * vec[I(x1,y1)]
						+ sx * ty * vec[I(x0,y1)];
	return val;
} 
