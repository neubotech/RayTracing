#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#pragma warning(disable: 4996)


#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include <time.h>
#include <math.h>

//#undef Success
#ifdef  _WIN32
#include "Eigen/Eigen/Core"
#include "Eigen/Eigen/Dense"
#else
#undef Success
#include <Eigen/Core>
#include <Eigen/Dense>
#endif //  _WIN32

#ifdef  _WIN32
#include "CImg/Cimg.h"
#else
#include "Cimg.h"
#endif // _WIN32


//****************************************************
// macros and inline functions
//****************************************************

#define PI 3.14159265  // Should be used from mathlib
#define EPS 1e-10       
#define DELETE_OBJECT(x)	if ((x)) { delete (x); (x) = NULL; }   // delete a class object
#define DELETE_ARRAY(x)		if ((x)) { delete [] (x); (x) = NULL; } // delete an array
#define FOR(i,n) for( int i=0; i<n; i++ )                           // for loop 
#define FOR_u(i, n) for (size_t i = 0; i < n; i++)                  // for loop 
#define SQUARE(x) ((x)*(x))
inline float sqr(float x) { return x*x; }
typedef unsigned char uchar; 
using namespace std;
using namespace Eigen;
using namespace cimg_library;


//****************************************************
// Basic Classes
//****************************************************
// class Vector{
// public:
// 	Vector(float x, float y, float z)
// public:
// };

// class Point:: public Vector3f{
// public:
// 	Point(float x, float y, float z)
// 	: m_coord(x, y, z){}

// public:
// 	Vector3f m_coord;
// };

// //all the m_coord are normalized to 1
// class Normal{
// public:
// 	Normal(float x, float y, float z)
// 	:m_coord(x,y,z){
// 		m_coord=m_coord/m_coord.norm();
// 	}
// public:
// 	Vector3f m_coord;
// 	// Vector4f m_hcoord;
// };

//***************************************************
//Ray
//*********************************************

class CRay{
public: 
	CRay(Vector3f pos, Vector3f dir, float t_min,float t_max)
	:m_pos(pos.x(), pos.y(), pos.z()),
	m_dir(dir.x(), dir.y(), dir.z()){
		m_dir=m_dir/m_dir.norm();

		m_t_min=t_min;
		m_t_max=t_max;

		// m_t=t_min;

		m_ray_start=pos+t_min*dir;
		m_ray_end=pos+t_max*dir;

		// m_ray=m_ray_start;
	}	

	Vector3f Ray_t(float t){
		if(t<m_t_min || t>m_t_max)
			return m_pos+t*m_dir;
		else
			cerr<<"t outside of range"<< m_t_min
				<<", "<< m_t_max <<endl;
			exit(-1);
	}

	bool hasPoint(Vector3f point){
		Vector3f dir_p=point-m_pos;

		//please check??
		return (dir_p-(dir_p).dot(m_dir)*m_dir).norm()==0;
	}


	bool hasPoint(Vector3f point, float threshold){
		Vector3f dir_p=point-m_pos;

		//please check??
		return (dir_p-(dir_p).dot(m_dir)*m_dir).norm()<=threshold;
	}

public:
	Vector3f m_pos, m_dir;
	Vector3f m_ray_start, m_ray_end;
	float m_t_min, m_t_max;

	// float m_t;
	// Vector3f m_ray;
};







//****************************************************
// color
//****************************************************

class CColor {
public: 
	CColor() { 
		m_rgb[0] = 0.0f; m_rgb[1] = 0.0f; m_rgb[2] = 0.0f; 
	}

	CColor(float _r, float _g, float _b) { 
		m_rgb[0] = _r;  m_rgb[1] = _g; m_rgb[2] = _b; 
	}
	
	CColor Add(CColor _c) {
		return CColor(m_rgb[0]+_c.m_rgb[0], m_rgb[1]+_c.m_rgb[1], m_rgb[2]+_c.m_rgb[2]);
	}

	void Print() {
		printf("(R, G, B) = (%3.3f, %3.3f, %3.3f)\n", m_rgb[0], m_rgb[1], m_rgb[2]);
	}

public: 
	float m_rgb[3];
};

//********************************************
//BRDF
//***************************************************
class CBRDF{
public:
	CColor ka, kd, ks, sp;
	CColor kr;
};


//****************************************************
// Light
//****************************************************
class CLight {
public: 
	enum ELightType {
		Point,
		Directional,
	};

	CLight() {}

	CLight(Vector3f _dir, CColor _color, ELightType _type) {
		m_dir = _dir; m_color = _color; m_type =  _type; 
	}

	void Print() {
		if (m_type == Point) 
			printf("Point Light ");
		else 
			printf("Directional Light ");
		printf("(x, y, z, R, G, B) = (%2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f)\n", 
			m_dir.x(), m_dir.y(), m_dir.z(), m_color.m_rgb[0], m_color.m_rgb[1], m_color.m_rgb[2]);
	}

public: 
	ELightType m_type;  // light type
	Vector3f m_dir;        //  location/position/direction
	CColor m_color;    // light color
};

//*************************************
//Camera (eye and screen and ray)
//******************************************
struct Pixel{
	Vector3f m_coord;
	CColor m_color;
}Pixel_t;


class CCamera {
public:
	CCamera(){}

	// : m_w(768), m_h(512), m_width(30.0), m_height(20.0),
	// 	m_center(0,0,0), m_normal(0,0,1){}

	CCamera(Vector3f _eye, int _w, int _h, Vector3f LL, Vector3f UL, Vector3f LR, Vector3f UR)
		:m_eye(_eye.x(), _eye.y(), _eye.z()),
		m_w(_w), m_h(_h), 
		m_width((LL-LR).norm()), m_height((UL-LL).norm()){
			m_center=((LL+LR)/2.0+(UL+UR)/2.0)/2.0;

			//make sure direction is right
			//normal pointing toward eye
			m_normal=(LL-LR).cross(LL-UL);
			m_normal=m_normal/m_normal.norm();
			
			m_pixel=new Pixel[m_w*m_h];

			float u_step=1.0/_w, v_step=1.0/_h;
			float u=u_step/2.0, v=v_step/2.0;

			for (int i=0; i<m_h; i++){
				for(int j=0; j< m_w; j++){
					m_pixel[i*m_w+j].m_coord=u*(v*LL+(1-v)*UL)+(1-u)*(v*LR+(1-v)*UR);
					m_pixel[i*m_w+j].m_color.m_rgb[0]=0.0;
					m_pixel[i*m_w+j].m_color.m_rgb[1]=0.0;
					m_pixel[i*m_w+j].m_color.m_rgb[2]=0.0;
					u+=u_step;
				}
				v+=v_step;
			}
		}

	void ColorPixel(int i, int j, CColor color){
		m_pixel[i*m_w+j].m_color.m_rgb[0]=color.m_rgb[0];
		m_pixel[i*m_w+j].m_color.m_rgb[1]=color.m_rgb[1];
		m_pixel[i*m_w+j].m_color.m_rgb[2]=color.m_rgb[2];
	}

	virtual void shoot(){}



	void DestroyCamera(){
		delete[] m_pixel;
	}

public:
	int m_w, m_h; // width and height  in pixels
	float m_width, m_height; //width and height in world units (cm)
	Vector3f m_center;
	Vector3f m_normal;
	Vector3f LL, UL, LR, UR;
	Vector3f m_eye;
	Pixel *m_pixel;

};


//****************************************************
// shape 
//****************************************************
class C3DModel {
public: 
	//virtual CShape() = 0;
	virtual Vector3f Norm(Vector3f _pos) = 0;
	virtual bool IsVisible(float _x, float _y) = 0; 
	virtual float Z(float _x, float _y) = 0;
	virtual Vector3f Center() = 0;
	virtual float Radius() = 0;
public: 

};

class CSphere : public C3DModel {
public: 
	CSphere() {}
	CSphere(Vector3f _center, float _radius) {
		m_center = _center; 
		m_radius = _radius;
	}

	bool IsVisible(float _x, float _y) {
		return SQUARE(_x-m_center.x()) + SQUARE(_y-m_center.y()) < SQUARE(m_radius);
	}

	float Z(float _x, float _y) {
		return sqrt(SQUARE(m_radius) - SQUARE(_x-m_center.x()) - SQUARE(_y-m_center.y()))+m_center.z(); 
	}

	float Radius() { return m_radius;}
	Vector3f Center() { return m_center; }
public: 
	Vector3f m_center; 
	float m_radius; 
};

// class Screen {
// public:
// 	Screen(){ screen.resize(256, 256); m_center[0]=128; m_center[1]=128;}
// 	Screen(int _w, int _h){
// 		m_w
// 	}
	
// public:
// 	int m_center[2];
// 	int m_w, m_h;

// 	pngwriter screen(1,1,0,);

// }



//****************************************************
// Global Variables
//****************************************************

CColor g_Ka, g_Kd, g_Ks;  // parameters 
vector<CLight> g_lights;
vector<C3DModel*> g_models; 
float g_v; 
bool g_isMultiple = false; 
bool g_isToon = false;  
bool g_isSave = true; 

Vector3f d; //anisotropic fibrics direction vector
bool g_isAniso = false;


//USING_PART_OF_NAMESPACE_EIGEN

int main(int argc, char *argv[]){
	Matrix3f m3;
	m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	Matrix4f m4=Matrix4f::Identity();

	Vector3f p_ray(0,0,0), d_ray(1,1,0), p_point(1.9, 1.9, 0);
	CRay ray(p_ray, d_ray, 1, 10);

	cout<< ray.m_ray_start<< "\n\n" << ray.m_ray_end
		<<"\n\n"<<p_point<<"\n\n"
		<<ray.hasPoint(p_point, .1)<<endl;
	// cout<<"m3\n" << m3 << "\nm4:\n"
	// 	<<m4 << "\nv4:\n" <<N.m_coord<<endl;


	// CImg<unsigned char> image("lena.pgm"), visu(500,400,1,3,0);
 //  const unsigned char red[] = { 255,0,0 }, green[] = { 0,255,0 }, blue[] = { 0,0,255 };
 //  image.blur(2.5);
 //  CImgDisplay main_disp(image,"Click a point"), draw_disp(visu,"Intensity profile");
 //  while (!main_disp.is_closed() && !draw_disp.is_closed()) {
 //    main_disp.wait();
 //    if (main_disp.button() && main_disp.mouse_y()>=0) {
 //      const int y = main_disp.mouse_y();
 //      visu.fill(0).draw_graph(image.get_crop(0,y,0,0,image.width()-1,y,0,0),red,1,1,0,255,0);
 //      visu.draw_graph(image.get_crop(0,y,0,1,image.width()-1,y,0,1),green,1,1,0,255,0);
 //      visu.draw_graph(image.get_crop(0,y,0,2,image.width()-1,y,0,2),blue,1,1,0,255,0).display(draw_disp);
 //      }
 //    }
  return 0;

}
