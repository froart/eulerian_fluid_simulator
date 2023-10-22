#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>
#include "../include/eulerian_fluid_dynamics.hpp"
#include <vt_framebuffer.hpp>

using namespace std;

#define U_FIELD 0
#define V_FIELD 1
#define S_FIELD 2

#define I(y,x) ((x)+(this->pParameters->nx)*(y)) // index of an image

struct Fields
{

  vector<float>  sm;   // smoke density
  vector<float>  sm_p; // previous smoke density
  vector<float>  vx;   // velocity on x-axis
  vector<float>  vy;   // velocity on y-axis
  vector<float>  vx_p; // previous velocity on x-axis
  vector<float>  vy_p; // previous velocity on y-axis
  vector<float>  pr;   // pressure
  vector<float>  oc;   // occupance of the cell by fluid
 
};

struct Parameters
{

  float cell_size;
  int   nx; // number of cells on x-axis
  int   ny; // number of cells on y-axis
  float dt; // time differential 
  int   niter; // number of iterations for backward finite difference
  float dens; // density 
  float ovrelax; // over relaxion coefficient

};

Gas::Gas( float cell_size, 
          int nx, 
          int ny, 
          float dt, 
          float dens, 
          int niter )
{

	int cell_num = nx * ny;

  this->pFields       = new Fields;
  this->pFields->sm   = move( vector<float>( cell_num, 0.0 ) );
	this->pFields->sm_p = move( vector<float>( cell_num, 0.0 ) );
	this->pFields->vx   = move( vector<float>( cell_num, 0.0 ) );
	this->pFields->vy   = move( vector<float>( cell_num, 0.0 ) );
	this->pFields->vx_p = move( vector<float>( cell_num, 0.0 ) );
	this->pFields->vy_p = move( vector<float>( cell_num, 0.0 ) );
	this->pFields->pr   = move( vector<float>( cell_num, 0.0 ) );
	this->pFields->oc   = move( vector<float>( cell_num, 1.0 ) ); // 1.0 for fluid / 0.0 for obstacle
  
  this->pParameters            = new Parameters;
  this->pParameters->cell_size = cell_size;
  this->pParameters->nx        = nx;
  this->pParameters->ny        = ny;
  this->pParameters->dt        = dt;
  this->pParameters->niter     = niter;
  this->pParameters->ovrelax   = 1.9;

	for( int j = 0; j < ny; ++j )
		for( int i = 0; i < nx; ++i )
    {
			if( i == 0 || j == 0 || j == ny - 1 || i == nx - 1 ) // borders
				this->pFields->oc[I( j, i )] == 0.0; // 0.0 for obstacle
		}
};

void Gas::addSmoke( int y, int x, float amount )
{

	this->pFields->sm[I( y, x )] = amount;
 
}

void Gas::addWind( int y, int x, float y_amount, float x_amount )
{

	this->pFields->vx[I( y, x )] = x_amount;
	this->pFields->vy[I( y, x )] = y_amount;

}

void Gas::evaluate()
{

	advect_velocity();
	advect_smoke();
	project_fluid();
	extrapolate();

}

void Gas::project_fluid()
{ // force imcompressibility

	for(int k = 0; k < this->pParameters->niter; k++)
		for(int j = 1; j < this->pParameters->ny-1; j++)
			for(int i = 1; i < this->pParameters->nx-1; i++)
      {
			  if(!this->pFields->oc[I( j, i )]) continue; // if evaluating at an obstacle
			  float sx0 = this->pFields->oc[I( j, i-1 )];
			  float sx1 = this->pFields->oc[I( j, i+1 )];
			  float sy0 = this->pFields->oc[I( j-1, i )];
			  float sy1 = this->pFields->oc[I( j+1, i )];
			  float s = sx0 + sx1 + sy0 + sy1;
			  if(!s) continue; // if nothing to compute here...
			  float div = + this->pFields->vx[I( j, i+1 )] 
                    - this->pFields->vx[I( j, i   )] 
                    + this->pFields->vy[I( j+1, i )] 
                    - this->pFields->vy[I( j, i   )];
			  // computing the pressure (not necessary for the simulation
			  // FIXME  should be minus???
			  float p = -(div / s) * this->pParameters->ovrelax;
			  // TODO where does this equation come from???
			  this->pFields->pr[I( j, i )] += p * (this->pParameters->dens * this->pParameters->cell_size) / this->pParameters->dt; 
			  // fix the fluid to be imcompressible
			  this->pFields->vx[I( j, i   )] -= sx0 * p;
			  this->pFields->vx[I( j, i+1 )] += sx1 * p;
			  this->pFields->vy[I( j, i   )] -= sy0 * p;
			  this->pFields->vy[I( j+1, i )] += sy1 * p;
			}

}

void Gas::extrapolate()
{ // enforce border conditions

	for(int j = 0; j < this->pParameters->ny; ++j)
  {
		this->pFields->vx[I( j, 0 )] = -this->pFields->vx[I( j, 1 )]; // left wall
		this->pFields->vx[I( j, this->pParameters->nx-1 )] = -this->pFields->vx[I( j, this->pParameters->nx-2 )]; // right wall
	}

	for(int i = 0; i < this->pParameters->nx; ++i)
  {
		this->pFields->vy[I( 0, i )] = -this->pFields->vy[I( 1, i )]; // bottom wall
		this->pFields->vy[I( this->pParameters->ny-1 , i )] = -this->pFields->vy[I( this->pParameters->ny-2, i )]; // top wall
	}

}

void Gas::advect_velocity()
{

	this->pFields->vx_p = this->pFields->vx;
	this->pFields->vy_p = this->pFields->vy;
	float h = this->pParameters->cell_size;
	float h2 = h * 0.5;
	for(int j = 1; j < this->pParameters->ny-1; j++)
		for(int i = 1; i < this->pParameters->nx-1; i++)
    {
			if( this->pFields->oc[I( j, i )]
          && this->pFields->oc[I( j-1, i )]
          && j < this->pParameters->ny-1 ) 
      { // x-component
				float x = static_cast<float>( i ) * h ;
				float y = static_cast<float>( j ) * h + h2; // x-component is situated at this point
				float u = this->pFields->vx[I( j, i )];	
				float v = ( + this->pFields->vy[I( j, i-1 )] 
                    + this->pFields->vy[I( j+1, i-1 )] 
                    + this->pFields->vy[I( j, i )] 
                    + this->pFields->vy[I( j+1, i)] ) * 0.25;	// y-component in this case is averaged by 4 values around
				x -= u * this->pParameters->dt;
				y -= v * this->pParameters->dt;
				u = sample_field( y, x, U_FIELD, this->pFields->vx );
				this->pFields->vx_p[I( j, i )] = u;
			}

			if( this->pFields->oc[I( j, i )] 
          && this->pFields->oc[I( j, i-1 )] 
          && i < this->pParameters->nx-1 ) 
      { // y-component
				float x = static_cast<float>( i ) * h + h2; // y-component is situated at this point
				float y = static_cast<float>( j ) * h;
				float u = ( + this->pFields->vx[I( j, i )]
                    + this->pFields->vx[I( j-1, i )]
                    + this->pFields->vx[I( j-1, i+1 )] 
                    + this->pFields->vx[I( j, i+1 )] ) * 0.25;	// x-component in this case is averaged by 4 values around
				float v = this->pFields->vy[I( j, i )];	
				x -= u * this->pParameters->dt;
				y -= v * this->pParameters->dt;
				v = sample_field(y, x, V_FIELD, this->pFields->vy );
				this->pFields->vy_p[I( j, i )] = v;
			}
		}
	this->pFields->vx = this->pFields->vx_p;
	this->pFields->vy = this->pFields->vy_p;

}

void Gas::advect_smoke()
{

	this->pFields->sm_p = this->pFields->sm;
	float h = this->pParameters->cell_size;
	float h2 = 0.5 * h;
	for(int j = 1; j < this->pParameters->ny-1; ++j)
		for(int i = 1; i < this->pParameters->nx-1; ++i)
			if(this->pFields->oc[I( j, i )] ) 
      {
				float u = ( this->pFields->vx[I( j, i )] + this->pFields->vx[I( j, i+1 )] ) * 0.5;	
				float v = ( this->pFields->vy[I( j, i )] + this->pFields->vy[I( j+1, i )] ) * 0.5;	
				float x = static_cast<float>( i ) * h + h2 - this->pParameters->dt * u;
				float y = static_cast<float>( j ) * h + h2 - this->pParameters->dt * v;
				this->pFields->sm_p[I( j, i )] = sample_field(y, x, S_FIELD, this->pFields->sm );
			}
	this->pFields->sm = this->pFields->sm_p;

}

float Gas::sample_field( float y_p, float x_p, int field, vector<float>& field_vec )
{

	float h  = this->pParameters->cell_size;
	float h1 = 1.0 / this->pParameters->cell_size;
	float h2 = this->pParameters->cell_size * 0.5;
	
	float x = fmax( fmin( x_p, this->pParameters->nx * h ), h );
	float y = fmax( fmin( y_p, this->pParameters->ny * h ), h );

	float dx = 0.0, dy = 0.0;
	switch( field )
  {
		case U_FIELD: dy = h2;          break;
		case V_FIELD: dx = h2;          break;
		case S_FIELD: dx = h2; dy = h2; break;
	}
	// TODO: what is it?
	float x0 = fmin( floor( ( x-dx ) * h1), this->pParameters->nx-1 ); 	
	float tx = ( ( x-dx ) - x0 * h ) * h1;
	float x1 = fmin( x0 + 1, this->pParameters->nx-1 );

	float y0 = fmin( floor( ( y-dy ) * h1 ), this->pParameters->ny-1 ); 	
	float ty = ( ( y-dy ) - y0 * h ) * h1;
	float y1 = fmin( y0 + 1, this->pParameters->ny-1 );

	float sx = 1.0 - tx;
	float sy = 1.0 - ty;

	float val = sx * sy * field_vec[I( y0, x0 )]
						+ tx * sy * field_vec[I( y0, x1 )]
						+ tx * ty * field_vec[I( y1, x1 )]
						+ sx * ty * field_vec[I( y1, x0 )];
	return val;

} 

vector<float>& Gas::getImage()
{

  return this->pFields->sm;

}

Gas::~Gas()
{

  delete this->pParameters;
  delete this->pFields;

}
