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

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>

#undef Success
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Cimg.h"

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
	Ray(Vector3f pos, Vector3f dir, float t_min,float t_max)
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


//**************************************************
//Screen
//******************************************************

class CScreen {
public:
	CScreen(): m_w(768), m_h(512), m_width(30.0), m_height(20.0),
		m_center(0,0,0), m_normal(0,0,1){}

public:
	int m_w, m_h; // width and height  in pixels
	float m_width, m_height; //width and height in world units (cm)
	Vector3f m_center;
	Vector3f m_normal;
};

//*************************************
//Camera (eye)
//******************************************
class CCamera{

public: 
	Vector3f m_eye;
	CScreen m_screen;
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

	CLight(CVec3 _dir, CColor _color, ELightType _type) {
		m_dir = _dir; m_color = _color; m_type =  _type; 
	}

	void Print() {
		if (m_type == Point) 
			printf("Point Light ");
		else 
			printf("Directional Light ");
		printf("(x, y, z, R, G, B) = (%2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f)\n", 
			m_dir.m_x, m_dir.m_y, m_dir.m_z, m_color.m_rgb[0], m_color.m_rgb[1], m_color.m_rgb[2]);
	}

public: 
	ELightType m_type;  // light type
	CVec3 m_dir;        //  location/position/direction
	CColor m_color;    // light color
};

//****************************************************
// shape 
//****************************************************
class C3DModel {
public: 
	//virtual CShape() = 0;
	virtual CVec3 Norm(CVec3 _pos) = 0;
	virtual bool IsVisible(float _x, float _y) = 0; 
	virtual float Z(float _x, float _y) = 0;
	virtual CVec3 Center() = 0;
	virtual float Radius() = 0;
public: 

};

class CSphere : public C3DModel {
public: 
	CSphere() {}
	CSphere(CVec3 _center, float _radius) {
		m_center = _center; 
		m_radius = _radius;
	}

	CVec3 Norm(CVec3 _pos) {
		return _pos.Subtract(m_center).UnitVec();
	}

	bool IsVisible(float _x, float _y) {
		return SQUARE(_x-m_center.m_x) + SQUARE(_y-m_center.m_y) < SQUARE(m_radius);
	}

	float Z(float _x, float _y) {
		return sqrt(SQUARE(m_radius) - SQUARE(_x-m_center.m_x) - SQUARE(_y-m_center.m_y))+m_center.m_z; 
	}

	float Radius() { return m_radius;}
	CVec3 Center() { return m_center; }
public: 
	CVec3 m_center; 
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

CScreen g_viewport;
CColor g_Ka, g_Kd, g_Ks;  // parameters 
vector<CLight> g_lights;
vector<C3DModel*> g_models; 
float g_v; 
bool g_isMultiple = false; 
bool g_isToon = false;  
bool g_isSave = true; 

CVec3 d; //anisotropic fibrics direction vector
bool g_isAniso = false;


//USING_PART_OF_NAMESPACE_EIGEN

int main(int argc, char *argv[]){
	Matrix3f m3;
	m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	Matrix4f m4=Matrix4f::Identity();

	Vector3f p_ray(0,0,0), d_ray(1,1,0);
	CRay ray(p_ray, d_ray, 1, 10);

	cout<< ray.m_ray_start<< "\n\n" << ray.m_ray_end
		<<"\n\n"<<ray.hasPoint(Vector3f(0.7,1.8,0), .1)<<endl;
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
