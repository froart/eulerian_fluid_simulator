#pragma once

#include <vector>

using namespace std;

class Gas {
  private:
    struct Fields* pFields;
    struct Parameters* pParameters;
	  void  project_fluid();
	  void  extrapolate();
	  void  advect_velocity();
	  void  advect_smoke();
	  float sample_field( float, float, int, vector<float>& );
	public:
	  Gas( vector<float>*, float, int, int, float, float, int );
	  void addSmoke(int, int, float);
	  void addWind(int, int, float, float);
	  void evaluate();
    ~Gas();
};

