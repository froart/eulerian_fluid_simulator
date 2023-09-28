#pragma once

#include <vector>

using namespace std;

class Fluid {
  private:
	  float m_cell_size;
	  int m_nx; // number of cells on x-axis
	  int m_ny; // number of cells on y-axis
	  float m_dt; // time differential 
	  int m_niter; // number of iterations for backward finite difference
	  float m_dens; // density 
	  float m_ovrelax; // over relaxion coefficient
	  vector<float> m_vx; // velocity on x-axis
	  vector<float> m_vy; // velocity on y-axis
	  vector<float> m_vx1; // previous velocity on x-axis
	  vector<float> m_vy1; // previous velocity on y-axis
	  vector<float> m_p; // pressure
	  vector<float> m_s; // occupance of the cell by fluid
	  vector<float> m_m; // smoke density
	  vector<float> m_m1; // previous smoke density
	  void project_fluid();
	  void extrapolate();
	  void advect_velocity();
	  void advect_smoke();
	  float sample_field(float, float, int, vector<float>&);
	public:
	  Fluid(float, int, int, float, float, int);
	  float* getImage();
	  void addSmoke(int, int, float);
	  void addWind(int, int, float, float);
	  void evaluate();
	  void clearImage();
};

