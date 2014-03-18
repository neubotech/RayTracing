#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
//#include "Timer.h"
//#include <openmp>
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
#include "Eigen/Eigen/Core"
#include "Eigen/Eigen/Dense"
#endif //  _WIN32

#ifdef  _WIN32
#include "CImg/Cimg.h"
#else
#include "Cimg.h"
#endif // _WIN32

using namespace std;
using namespace Eigen;
using namespace cimg_library;


//****************************************************
// macros and inline functions
//****************************************************
#define EPS 1e-5f
#define EPS_SELF 1e-6f
#define NUM_DEPTH 2
//#define NUM_LIGHTS 4
//#ifdef _WIN32
//#define INF 1e10f
//#endif // _WIN32

#define DELETE_OBJECT(x)	if ((x)) { delete (x); (x) = NULL; }   // delete a class object
#define DELETE_ARRAY(x)		if ((x)) { delete [] (x); (x) = NULL; } // delete an array
#define FOR(i,n) for( int i=0; i<n; i++ )                           // for loop 
#define FOR_u(i, n) for (size_t i = 0; i < n; i++)                  // for loop 
#define SQUARE(x) ((x)*(x))
#define INF (float) 1e10
#define PI (float) 3.1415926

inline float sqr(float x) { return x*x; }
typedef unsigned char uchar; 
typedef Vector3f V3f; 
#define PRINT
//***************************************************
//Ray
//*********************************************

class CRay{
public: 
	CRay() {}
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
		if(t >= m_t_min && t <= m_t_max)
			return m_pos+t*m_dir;
		else {
			cerr<<"t outside of range"<< m_t_min
				<<", "<< m_t_max <<endl;
			exit(-1);
		}
	}

	/*void Print() {
		cout << m_dir;
	}*/

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

	void Print(){
		#ifdef PRINT
		cout<<"[Ray] pos:\n"<<m_pos<<"\ndir:\n"<<m_dir
			<< "\ntmin: " << m_t_min << " tmax: " << m_t_max << endl; 
		#endif // PRINT

	
			//<<" start: "<<m_ray_start<<" end: "<< m_ray_start<<endl;
	}

public:
	Vector3f m_pos, m_dir;
	Vector3f m_ray_start, m_ray_end;
	float m_t_min, m_t_max;

	// float m_t;
	// Vector3f m_ray;
};



//****************************************************
// LocalGeo
//****************************************************
//Notes: Store the local geometry at the intersection point. May need to store
//		other quantities (eg. texture coordinate) in a more complicated raytracer.
class CLocalGeo {
public: 
	Vector3f m_pos;   // point 
	Vector3f m_n;     // normal (unit vector)
};

//****************************************************
// shape
//****************************************************
void TestShape() {

}

class CShape {
public:
	virtual bool intersect(CRay& _ray, float* _thit, CLocalGeo* _local) = 0; 
	virtual bool intersectP(CRay& _ray) = 0; 
};

class CSphere : public CShape {
public: 
	CSphere(V3f _c, float _r) {
		m_c = _c;  m_r = _r; m_r2 = m_r * m_r; 
	}

	bool intersect(CRay& _ray, float* _thit, CLocalGeo* _local)  {
		float t; 
		V3f oc = _ray.m_pos - m_c; 
		float a = _ray.m_dir.dot(_ray.m_dir); // dir should be unit
		float b = _ray.m_dir.dot(oc);
		float c = oc.dot(oc) - m_r2; 
		float delta = b * b  - a * c; 
		if (delta < 0)   // no solution 
			return false; 
		else if (delta > -EPS && delta < EPS) {  // one solution
			t = - b / a;
			if (t > _ray.m_t_max || t < _ray.m_t_min)  // out of range
				return false; 
		} else {   // two solutions 
			float deltasqrt = sqrt(delta);
			float t1 = (- b - deltasqrt) / a;
			float t2 = (- b + deltasqrt) / a;
			t = _ray.m_t_max + 1; 
			if (t1 <= _ray.m_t_max && t1 >= _ray.m_t_min)
				t = min(t, t1);

			if (t2 <= _ray.m_t_max && t2 >= _ray.m_t_min)
				t = min(t, t2);
			if (t > _ray.m_t_max)   // both out of range
				return false; 
		}

		// pass t, compute CLocalGeo
		*_thit = t; 
		_local->m_pos = _ray.Ray_t(t);
		_local->m_n = _local->m_pos - m_c; 
		_local->m_n = _local->m_n / _local->m_n.norm(); 
		return true; 
	} 
	
	bool intersectP(CRay& _ray) {
		V3f oc = _ray.m_pos - m_c; 
		float a = _ray.m_dir.dot(_ray.m_dir); // dir should be unit
		float b = _ray.m_dir.dot(oc);
		float c = oc.dot(oc) - m_r2; 
		float delta = b * b  - a * c; 
		if (delta < 0)   // no solution 
			return false; 
		else if (delta > -EPS && delta < EPS) {  // one solution
			float t = - b / a;
			if (t > _ray.m_t_max || t < _ray.m_t_min)  // out of range
				return false; 
		} else {   // two solutions 
			float deltasqrt = sqrt(delta);
			float t1 = (- b - deltasqrt) / a;
			float t2 = (- b + deltasqrt) / a;
			if (t1 > _ray.m_t_max || t1 < _ray.m_t_min)
				return false; 

			if (t2 > _ray.m_t_max || t2 < _ray.m_t_min)
				return false; 
		}

		return true; 
	}

private: 
  	V3f m_c;     // center
	float m_r;  // radius
	float m_r2; // radius * radius;
};


class CTriangle : public CShape {
private: 
	float det(V3f _V0, V3f _V1, V3f _V2) {
		float d = _V0.x() * _V1.y() * _V2.z() 
			+ _V1.x() * _V2.y() * _V0.z()
			+ _V2.x() * _V0.y() * _V1.z()
			- _V2.x() * _V1.y() * _V0.z()
			- _V1.x() * _V0.y() * _V2.z()
			- _V0.x() * _V2.y() * _V1.z();
		return d; 
	}

public: 
	CTriangle(V3f _V0, V3f _V1, V3f _V2) {
		m_V0 = _V0; m_V1 = _V1; m_V2 = _V2; 
	}

	//bool intersectP(CRay& _ray) {
	//	V3f E1 = m_V0 - m_V1; 
	//	V3f E2 = m_V0 - m_V2; 
	//	V3f S = m_V0 - _ray.m_pos;
	//	float d = det(_ray.m_dir, E1, E2);
	//	if (d < EPS && d > -EPS)
	//		return false; 
	//	float d1 = det(S, E1, E2);
	//	float t = d1 / d; 
	//	if (t < _ray.m_t_min || t > _ray.m_t_max)
	//		return false; 
	//	float d2 = det(_ray.m_dir, S, E2);
	//	float beta = d2 / d; 
	//	if (beta < -EPS || beta - 1 > EPS)
	//		return false; 
	//	float d3 = det(_ray.m_dir, E1, S);
	//	float gamma = d3 / d; 
	//	if (gamma < -EPS || gamma - 1 > EPS || beta + gamma - 1 > EPS)
	//		return false; 
	//	return true; 
	//}

	//bool intersect(CRay& _ray, float* _thit, CLocalGeo* _local)  {
	//	V3f E1 = m_V0 - m_V1; 
	//	V3f E2 = m_V0 - m_V2; 
	//	V3f S = m_V0 - _ray.m_pos;
	//	float d = det(_ray.m_dir, E1, E2);
	//	if (d < EPS && d > -EPS)
	//	 	return false; 
	//	float d1 = det(S, E1, E2);
	//	float t = d1 / d; 
	//	if (t < _ray.m_t_min || t > _ray.m_t_max)
	//	 	return false; 
	//	float d2 = det(_ray.m_dir, S, E2);
	//	float beta = d2 / d; 
	//	if (beta < -EPS || beta  - 1 > EPS)
	//	 	return false; 
	//	float d3 = det(_ray.m_dir, E1, S);
	//	float gamma = d3 / d; 
	//	if (gamma < -EPS || gamma - 1 > EPS || beta + gamma - 1 > EPS)
	//	 	return false; 
	//	*_thit = t;
	//	_local->m_pos = _ray.Ray_t(t);
	//	_local->m_n = (m_V0-m_V1).cross(m_V0-m_V2);
	//	_local->m_n = _local->m_n / _local->m_n.norm();
	//	return true; 
	//}
	bool intersect(CRay& _ray, float* _thit, CLocalGeo* _local)  {
		V3f A = m_V1 - m_V0; 
		V3f B = m_V2 - m_V0; 
		// cout<<A<<" , "<<B<<endl;
		V3f N=A.cross(B);
		float notParallel=N.dot(_ray.m_dir);
		if(notParallel==0){
			//cout<<"[triangle] ray parallel"<<endl;
			return false; //triangle is parallel
		}
		float d=N.dot(m_V0);
		float t=(N.dot(_ray.m_pos)+d)/notParallel;


		if(t<_ray.m_t_min || t > _ray.m_t_max){
			//cout<< "[triangle] ray is behind"<<endl;
			return false; //triangle is behind
		}

		V3f P=_ray.m_pos + t*_ray.m_dir;
		V3f C;

		V3f edge0 = m_V1 - m_V0;
		V3f VP0 = P - m_V0;
		C=edge0.cross(VP0);
		if(N.dot(C)<0){
			//cout<<"on the right of edge0"<<endl;
			return false; //P is on the right side
		}

		V3f edge1= m_V2 - m_V1;
		V3f VP1 = P - m_V1;
		C=edge1.cross(VP1);
		if(N.dot(C)<0){
			//cout<<"on the right of edge1"<<endl;
			return false; //P is on the right side
		}

		V3f edge2= m_V0 - m_V2;
		V3f VP2 = P - m_V2;
		C=edge2.cross(VP2);
		if(N.dot(C)<0){
			//cout<<"on the right of edge2"<<endl;
			return false; //P is on the right side
		}
		//cout<<"is inside"<<endl;

		*_thit = t;
		_local->m_pos = _ray.Ray_t(t);
		_local->m_n = (m_V0-m_V1).cross(m_V0-m_V2);
		_local->m_n = _local->m_n / _local->m_n.norm();
		return true; 
	} 

	// bool intersect(CRay& _ray, float* _thit, CLocalGeo* _local)  {
	// 	V3f E1 = m_V0 - m_V1; 
	// 	V3f E2 = m_V0 - m_V2; 
	// 	V3f S = m_V0 - _ray.m_pos;
	// 	float d = det(_ray.m_dir, E1, E2);
	// 	if (d < EPS && d > -EPS)
	// 		return false; 
	// 	float d1 = det(S, E1, E2);
	// 	float t = d1 / d; 
	// 	if (t < _ray.m_t_min || t > _ray.m_t_max)
	// 		return false; 
	// 	float d2 = det(_ray.m_dir, S, E2);
	// 	float beta = d2 / d; 
	// 	if (beta < 0 || beta > 1)
	// 		return false; 
	// 	float d3 = det(_ray.m_dir, E1, S);
	// 	float gamma = d3 / d; 
	// 	if (gamma < 0 || gamma > 1 || beta + gamma > 1)
	// 		return false; 

	// 	*_thit = t;
	// 	_local->m_pos = _ray.Ray_t(t);
	// 	_local->m_n = m_V0.cross(m_V1);
	// 	_local->m_n = _local->m_n / _local->m_n.norm();
	// 	return true; 
	// } 

	bool intersectP(CRay& _ray) {
		V3f A = m_V1 - m_V0; 
		V3f B = m_V2 - m_V0; 
		// cout<<A<<" , "<<B<<endl;
		V3f N=A.cross(B);
		float notParallel=N.dot(_ray.m_dir);
		if(notParallel==0){
			//cout<<"[triangle] ray parallel"<<endl;
			return false; //triangle is parallel
		}
		float d=N.dot(m_V0);
		float t=(N.dot(_ray.m_pos)+d)/notParallel;


		if(t<_ray.m_t_min || t > _ray.m_t_max){
			//cout<< "[triangle] ray is behind"<<endl;
			return false; //triangle is behind
		}

		V3f P=_ray.m_pos + t*_ray.m_dir;
		V3f C;

		V3f edge0 = m_V1 - m_V0;
		V3f VP0 = P - m_V0;
		C=edge0.cross(VP0);
		if(N.dot(C)<0){
			//cout<<"on the right of edge0"<<endl;
			return false; //P is on the right side
		}

		V3f edge1= m_V2 - m_V1;
		V3f VP1 = P - m_V1;
		C=edge1.cross(VP1);
		if(N.dot(C)<0){
			//cout<<"on the right of edge1"<<endl;
			return false; //P is on the right side
		}

		V3f edge2= m_V0 - m_V2;
		V3f VP2 = P - m_V2;
		C=edge2.cross(VP2);
		if(N.dot(C)<0){
			//cout<<"on the right of edge2"<<endl;
			return false; //P is on the right side
		}

		//*_thit = t;
		//_local->m_pos = _ray.Ray_t(t);
		//_local->m_n = m_V0.cross(m_V1);
		//_local->m_n = _local->m_n / _local->m_n.norm();
		return true; 
	}
	//	//V3f E1 = m_V0 - m_V1; 
	//	//V3f E2 = m_V0 - m_V2; 
	//	//V3f S = m_V0 - _ray.m_pos;
	//	//float d = det(_ray.m_dir, E1, E2);
	//	//if (d < EPS && d > -EPS)   // zero 
	//	//	return false; 
	//	//float d1 = det(S, E1, E2);
	//	//float t = d1 / d; 
	//	//if (t < _ray.m_t_min || t > _ray.m_t_max)  // check t range 
	//	//	return false; 
	//	//float d2 = det(_ray.m_dir, S, E2);
	//	//float beta = d2 / d; 
	//	//if (beta < 0 || beta > 1)                 // check beta
	//	//	return false; 
	//	//float d3 = det(_ray.m_dir, E1, S);
	//	//float gamma = d3 / d;                       
	//	//if (gamma < 0 || gamma > 1 || beta + gamma > 1)   // check gamma
	//	//	return false; 
	//	//
	//	//return true; 
	//}

private: 
	V3f m_V0, m_V1, m_V2;     // 3D points
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
//BRDF and Material
//***************************************************
class CBRDF{
public: 
	CBRDF() {
		ka = CColor(0.0f, 0.0f, 0.0f);
		kd = CColor(0.0f, 0.0f, 0.0f);
		ks = CColor(0.0f, 0.0f, 0.0f);
		p = 0.0f; 
	}
public:
	CColor ka, kd, ks;
	CColor kr;
	float p; 
};

//Notes: Class for storing material. For this example, it just returns a constant
//	material regardless of what local is. Later on, when we want to support
//	texture mapping, this need to be modified.
class CMaterial {
public: 
	CMaterial(CBRDF _brdf) {
		m_constantBRDF = _brdf; 
	}
	CMaterial() {
		// set constant BRDF 
	}

	void getBRDF(CLocalGeo& _local, CBRDF* _brdf) {
		*_brdf = m_constantBRDF;
	}

private: 
	CBRDF m_constantBRDF;
};


//****************************************************
// primitive
//****************************************************
class CPrimitive; 

class CIntersection {
public: 
	CLocalGeo m_localGeo;
	CPrimitive* m_prim; 
};

class CMatrix
{
public:
	CMatrix(){
		m_translate<< 1, 0, 0, 0,
					0, 1, 0, 0,
					0, 0, 1, 0,
					0, 0, 0, 1;

		m_rotate = m_translate;
		m_scale = m_translate;

	}
	// Matrix(Matrix4f _mat4){
	// 	if(_mat4(15)!=0){
	// 		_mat4/=_mat4(15);
	// 	}
	// 	else{
	// 		cerr<<"matrix4f needs to be homogenious coordinate, mat[15] can't be 0";
	// 		exit(-1);
	// 	}

	// 	FOR (i, 16){ m_matrix4(i)=_mat4(i); }
	// 	m_matrix3<< m_matrix4(0), m_matrix4(1), m_matrix4(2),
	// 				m_matrix4(4), m_matrix4(5), m_matrix4(6),
	// 				m_matrix4(8), m_matrix4(9), m_matrix4(10);
	// }

	// Matrix(Matrix3f _mat3){
	// 	// FOR(i, 16){ m_matrix3(i)=_mat3(i); }

	// 	m_matrix4 << _mat3(0), _mat3(1), _mat3(2), 0.0f,
	// 				_mat3(3), _mat3(4), _mat3(5), 0.0f,
	// 				_mat3(6), _mat3(7), _mat3(8), 0.0f,
	// 					0.0f,	0.0f, 		0.0f, 1.0f;
	// }
	// void sync_3to4(){
	// 	m_matrix4 << m_matrix3(0), m_matrix3(1), m_matrix3(2), 0.0f,
	// 				m_matrix3(3), m_matrix3(4), m_matrix3(5), 0.0f,
	// 				m_matrix3(6), m_matrix3(7), m_matrix3(8), 0.0f,
	// 					0.0f,	0.0f, 		0.0f, 1.0f;
	// }

	// void sync_4to3(){

	// 	if(m_matrix4(15)!=0){
	// 		m_matrix4 /= m_matrix4(15); 
	// 		m_matrix3<< m_matrix4(0), m_matrix4(1), m_matrix4(2),
	// 					m_matrix4(4), m_matrix4(5), m_matrix4(6),
	// 					m_matrix4(8), m_matrix4(9), m_matrix4(10);
	// 			}
	// 	else{
	// 		cerr<<"matrix4f needs to be homogenious coordinate, mat[15] can't be 0";
	// 		exit(-1);
	// 	}
	// }

	//around axis, and counter clockwise (radius)
	void rotate_axis(Vector3f axis, float angle){
		angle*=(PI/180);
		axis/=axis.norm();
		Matrix3f r_cross;
		r_cross << 0, -axis(2), axis(1),
					axis(2), 0, -axis(0),
					-axis(1), axis(0), 0;

		Matrix3f temp= r_cross*r_cross.transpose()+sin(angle)*r_cross-cos(angle)*r_cross*r_cross;
		m_rotate<<temp(0), temp(1), temp(2), 0.0f,
				temp(4), temp(5), temp(6), 0.0f,
				temp(7), temp(8), temp(9), 0.0f,
				0.0f, 0.0f, 0.0f, 1.0f;

		// sync_3to4();
	}


	 // Matrix rotate_quat(){}
	//translate by distance in its 3 coordiate values
	void translate(Vector3f distance){
		m_translate << 1.0f, 0.0f, 0.0f, distance(0),
					0.0f, 1.0f, 0.0f, distance(1),
					0.0f, 0.0f, 1.0f, distance(2),
					0.0f, 0.0f, 0.0f, 1.0f;


		// m_matrix4 = temp * m_matrix4;
		// m_matrix4 = m_matrix4/m_matrix4(15);
		//m_matrix4 = temp;
		// sync_4to3();
	}

	//scale base on the x, y, z of vector 
	void scale(Vector3f factor){
		// Matrix3f temp;
		m_scale << factor(0), 	0.0f	, 0.0f	, 0.0f,
					0.0f, factor(1)	, 0.0f	,0.0f,
					0.0f, 0.0f		, factor(2), 0.0f,
					0.0f, 0.0f, 0.0f, 1;
		// m_matrix3 = temp*m_matrix3;
		// sync_3to4();
	}

	void inverse(){
		m_rotate=m_rotate.inverse();
		m_translate = m_translate.inverse();
		m_translate = m_translate/m_scale(15);
		m_scale = m_scale.inverse();
	}


public:
	// Matrix4f m_matrix4;
	// Matrix3f m_matrix3;

	Matrix4f m_translate, m_rotate, m_scale;

};

class CTransformation {

public: 
	CTransformation()
	{

		CMatrix m_Trans;
	}

	CTransformation( CMatrix _Trans ){
		m_Trans.m_rotate=_Trans.m_rotate;
		m_Trans.m_translate=_Trans.m_translate;
		m_Trans.m_scale=_Trans.m_scale;
	}


	V3f Dir(V3f _v) {
		Vector4f _v4;
		_v4 << _v(0), _v(1), _v(2), 0.0f;
		_v4 = m_Trans.m_rotate*_v4;
		_v << _v4(0), _v4(1), _v4(2);
		return _v; 
	}
	V3f Point(V3f _p) { 
		Vector4f _p4;
		_p4 << _p(0), _p(1), _p(2), 1.0f;
		_p4 = m_Trans.m_translate*m_Trans.m_rotate*m_Trans.m_scale*_p4;
		_p4 = _p4 / _p4(3);
		_p << _p4(0), _p4(1), _p4(2);
		return _p; 
	}

	V3f Normal(V3f _n) { 
		Vector4f _n4;
		_n4 << _n(0), _n(1), _n(2), 0.0f;
		_n4 = (m_Trans.m_rotate*m_Trans.m_scale).inverse().transpose()*_n4;
		_n4 = _n4 / _n4(3);
		_n << _n4(0), _n4(1), _n4(2);

		return _n; 
	}

	CRay Ray(CRay _ray) { 
		CRay n_ray;

		n_ray.m_pos = Point(_ray.m_pos);
		n_ray.m_dir = Dir(_ray.m_dir);

		n_ray.m_ray_start=Point(_ray.m_ray_start);
		n_ray.m_ray_end=Point(_ray.m_ray_end);

		n_ray.m_t_min=(n_ray.m_ray_start-n_ray.m_pos).norm();
		n_ray.m_t_max=(n_ray.m_ray_end-n_ray.m_pos).norm();

		return n_ray; 
	}
	CLocalGeo LocalGeo(CLocalGeo _localGeo) { 
		_localGeo.m_pos = Point(_localGeo.m_pos);
		_localGeo.m_n=Normal(_localGeo.m_n);
		return _localGeo; 
	}
public:
	CMatrix m_Trans;
};



class CPrimitive {
public: 
	virtual bool intersect(CRay& _ray, float* _thit, CIntersection* _in) = 0;
	virtual bool intersectP(CRay& _ray) = 0;
	virtual void getBRDF(CLocalGeo& local, CBRDF* brdf) = 0;
};

class CGeometricPrimitive : public CPrimitive {
public: 
	bool intersect(CRay& _ray, float* _thit, CIntersection* _in) {
		CRay oray = m_worldToObj.Ray(_ray);
		CLocalGeo olocal; 
		if (!m_shape->intersect(oray, _thit, &olocal)) 
			return false;
		_in->m_prim = this;
		_in->m_localGeo = m_objToWorld.LocalGeo(olocal);
		return true; 
	}

	bool intersectP(CRay& _ray) {
		CRay oray = m_worldToObj.Ray(_ray);
		return m_shape->intersectP(oray); 
	}

	void getBRDF(CLocalGeo& local, CBRDF* brdf) {
		m_mat->getBRDF(local, brdf);
	}

public: 
	CTransformation m_objToWorld; 
	CTransformation m_worldToObj; 
	CShape* m_shape; 
	CMaterial* m_mat; 
};


class CAggregatePrimitive : public CPrimitive {
public: 
	CAggregatePrimitive(vector<CPrimitive*> _primVec) {
		m_primVec = _primVec; 
		m_numPrims = (int)m_primVec.size();
	}

	bool intersect(CRay& _ray, float* _thit, CIntersection* _in) {
		float min_t = _ray.m_t_max + 1; 
		//CIntersection* min_in = NULL; 
		bool flag = false; 
		FOR (i, m_numPrims) {
			float t = 1e10f; 
			CIntersection in; // = new CIntersection(); 
			if (m_primVec[i]->intersect(_ray, &t, &in)) {
				flag = true; 
				if (t <= min_t) {
					min_t = t; 
					//DELETE_OBJECT(min_in);  //BUG?
					//min_in = _in;
					//DELETE_OBJECT(_in);
					//_in = in;
					//_in->m_localGeo = in->m_localGeo;
					//_in->m_prim = in->m_prim;
					//DELETE_OBJECT(in);
					*_in = in; 
				}
			}
		}

		return flag; 
	}

	bool intersectP(CRay& _ray) {
		bool flag = false; 
		FOR (i, m_numPrims)
			flag = flag || m_primVec[i]->intersectP(_ray);
		return flag; 
	}
	
	// This should never get called, because in->primitive will
	// never be an aggregate primitive
	void getBRDF(CLocalGeo& _local, CBRDF* _brdf) {
		printf("getBRDF should never get called\n");
		exit(-1);
	}

private: 
	int m_numPrims; 
	vector<CPrimitive*> m_primVec; 

};

// TO-DO 
class CBSPPrimitive : public CPrimitive {

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

	CLight(V3f _dir, CColor _color, ELightType _type) {
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

	// normalize light dir vector chech if direction is right????
	void generateLightRay(CLocalGeo& local, CRay* lray, CColor* lcolor) {
		*lcolor = m_color;
		if (m_type == Point) {
			lray->m_pos = local.m_pos + EPS_SELF * local.m_n; 
			lray->m_dir = m_dir - local.m_pos;
			lray->m_t_min = 0.0f; //1e-6;
			lray->m_t_max = lray->m_dir.norm();
			lray->m_dir = lray->m_dir/lray->m_dir.norm();
			
			lray->m_ray_start = lray->m_pos;
			lray->m_ray_end = m_dir;
		} 

		if (m_type == Directional) {
			lray->m_pos = local.m_pos + EPS_SELF * local.m_n; 
			lray->m_dir = m_dir / m_dir.norm();

			lray->m_t_min = 0.0f; 
			lray->m_t_max = INF; 
			lray->m_ray_start = lray->m_pos; 
			lray->m_ray_end = V3f(INF, INF, INF); //local.m_pos;

			// lray->m_t_min=0;
			// lray->m_t_max=INF;

			// lray->m_ray_start=local.m_pos;
			// lray->m_ray_end=INF;
		}
	}

public: 
	ELightType m_type;  // light type
	Vector3f m_dir;        //  location/position/direction
	CColor m_color;    // light color
};


class CRayTracer {

private:  
	CPrimitive* m_scene;     // primitives
	vector<CLight*> m_lights;   // lights

private: 
	CColor renderAmbientTerm(CRay& _ray, CColor& _rayColor, CColor& _Ka) {
		CColor c;
		FOR (k, 3) 
			c.m_rgb[k] = _Ka.m_rgb[k] * _rayColor.m_rgb[k];
		return c; 
	}

	CColor renderDiffuseTerm(CRay& _ray, CColor& _rayColor, CColor& _Kd, V3f& _n) {
		CColor c;
		float dot = _n.dot(_ray.m_dir);
		/*if (dot > 0.0) {
			_light.m_dir.UnitVec().Print();
			_n.Print();
		}*/
		dot = max(dot, 0.0f);
		//printf("dot = %2.2f\n", dot);
		FOR (k, 3)
			c.m_rgb[k] = _Kd.m_rgb[k] * dot * _rayColor.m_rgb[k];

		return c; 
	}

	CColor renderSpecularTerm(CRay& _ray, CColor& _rayColor, CColor& _Ks, V3f& _n, V3f& _v, float _p) {
		CColor c;
		//_ray.m_dir = _ray.m_dir;
		V3f dir = _ray.m_dir; 
		V3f r = (_n * 2 *dir.dot(_n))- dir;
		_v = -_v / _v.norm();

		float dot = pow(max(r.dot(_v),0.0f), _p);
		FOR (k, 3) 
			c.m_rgb[k] = _Ks.m_rgb[k] * dot * _rayColor.m_rgb[k];

		return c; 
	}

	CColor shading(CLocalGeo& _local, CBRDF& _brdf, CRay& _ray, V3f& _v, CColor& _rayColor) {
		CColor c_all(0.0f, 0.0f, 0.0f);

	/*	CColor c_ambient = renderAmbientTerm(_ray, _rayColor, _brdf.ka); 
		c_all = c_all.Add(c_ambient);*/

		CColor c_diffuse = renderDiffuseTerm(_ray, _rayColor, _brdf.kd, _local.m_n);
		c_all = c_all.Add(c_diffuse);

		CColor c_specular = renderSpecularTerm(_ray, _rayColor, _brdf.ks, _local.m_n, _v, _brdf.p);
		c_all = c_all.Add(c_specular);
		//c_all.Print(); 
		return c_all; 
	}

	CRay createReflectRay(CLocalGeo& _local, CRay& _ray) {
		V3f n = _local.m_n;
		V3f r = (n * 2 *_ray.m_dir.dot(n))- _ray.m_dir;
		//V3f pos = _local.m_pos + n * EPS_SELF; 
		V3f pos = _local.m_pos; 
		return CRay(pos, r, 1e-4f, INF);
	}

public: 
	CRayTracer() {}

	void Setup(CPrimitive* _scene, vector<CLight*> _lights) { // set up scene and lights
		m_scene = _scene; 
		m_lights = _lights;
	}

	bool trace(CRay& ray, int depth, CColor* color, CBRDF& brdf) {
		if (depth > NUM_DEPTH) {
			// Make the color black and return
			color->m_rgb[0] = 0.0; 
			color->m_rgb[1] = 0.0; 
			color->m_rgb[2] = 0.0; 
			return false; 
		}

		float thit; 
		CIntersection in;
		if (!m_scene->intersect(ray, &thit, &in)) {
			color->m_rgb[0] = 0.0; 
			color->m_rgb[1] = 0.0; 
			color->m_rgb[2] = 0.0; 
			return false;
		}

		//CBRDF brdf; 
		in.m_prim->getBRDF(in.m_localGeo, &brdf);

		FOR_u (i, m_lights.size()) {
			CRay lray; 
			CColor lcolor;
			//cout << in.m_localGeo.m_pos << endl; 
			//cout << in.m_localGeo.m_n << endl; 
			m_lights[i]->generateLightRay(in.m_localGeo, &lray, &lcolor);
			//lray.Print();
			if (true || !m_scene->intersectP(lray)) { //If not, do shading calculation for this light source
				// bug
				//shading(CLocalGeo& _local, CBRDF& _brdf, CRay& _ray, V3f& _v, CColor& _rayColor) {
				CColor c = shading(in.m_localGeo, brdf, lray, ray.m_dir, lcolor);
				*color = color->Add(c);
			}
		}

		// Handle mirror reflection
		if (brdf.kr.m_rgb[0] > 0 
			|| brdf.kr.m_rgb[1] > 0 
			|| brdf.kr.m_rgb[2] > 0) {
			CRay reflectRay = createReflectRay(in.m_localGeo, ray);
			//Make a recursive call to trace the reflected ray
			CColor tempColor; 
			CBRDF tempBRDF; 
			trace(reflectRay, depth+1, &tempColor,tempBRDF);
			color->m_rgb[0] += tempColor.m_rgb[0] * brdf.kr.m_rgb[0];
			color->m_rgb[1] += tempColor.m_rgb[1] * brdf.kr.m_rgb[1];
			color->m_rgb[2] += tempColor.m_rgb[2] * brdf.kr.m_rgb[2];
		}
		return true; 
	}
};

// global variables
vector<CLight*> g_lights; 
CPrimitive* g_scene; 
int m_width; 
int m_height; 
string g_fname; 



//*************************************
//Camera (eye and screen and ray)
//******************************************
struct Pixel{
	Vector3f m_coord;
	CColor m_color;
}Pixel_t;


class CCamera {
public:
	enum SamplerType { OverS, JitterS};
	CCamera(){}

	// : m_w(768), m_h(512), m_width(30.0), m_height(20.0),
	// 	m_center(0,0,0), m_normal(0,0,1){}
	void SetupRayTracer(CRayTracer* _tracer ) { m_rayTracer = _tracer; }
	CCamera(Vector3f _eye, int _w, int _h, Vector3f _LL, Vector3f _UL, Vector3f _LR, Vector3f _UR)
		:m_eye(_eye.x(), _eye.y(), _eye.z()),
		m_w(_w), m_h(_h), 
		m_width((LL-LR).norm()), m_height((UL-LL).norm()),
		LL(_LL.x(), _LL.y(), _LL.z()), UL(_UL.x(), _UL.y(), _UL.z()),
		LR(_LR.x(), _LR.y(), _LR.z()), UR(_UR.x(), _UR.y(), _UR.z()){
			m_center=((LL+LR)/2.0f+(UL+UR)/2.0f)/2.0f;

			//make sure direction is right
			//normal pointing toward eye
			m_normal=(LL-LR).cross(LL-UL);
			m_normal=m_normal/m_normal.norm();

			// m_sample=new Vector3f[m_w*m_h*m_over_sample_ratio*m_over_sample_ratio];
			m_pixel=new Pixel[m_w*m_h];

			//////////////////////////////////////////
			//pixel coordinates and color initialization
			//////////////////////////////////////////


			float u_step=1.0f/m_w;
			float v_step=1.0f/m_h;

			float u=u_step/2.0f;
			float v=v_step/2.0f;

			for (int i=0; i<m_h; i++){
				for(int j=0; j< m_w; j++){
					m_pixel[i*m_w+j].m_coord=u*(v*UR+(1-v)*LR)+(1-u)*(v*UL+(1-v)*LL);
					// m_pixel[i*m_w+j].m_color.m_rgb[0]=0.0;
					// m_pixel[i*m_w+j].m_color.m_rgb[1]=0.0;
					// m_pixel[i*m_w+j].m_color.m_rgb[2]=0.0;

					// cout<<m_pixel[i*m_w+j].m_coord<<endl;
					u+=u_step;
				}
				u=u_step/2.0f;
				v+=v_step;
			}
	}

	void ColorPixel(int i, int j, CColor color){
		float mag = 1.0f; 
		//color.Print(); 
		m_pixel[i*m_w+j].m_color.m_rgb[0]=color.m_rgb[0]*mag;
		m_pixel[i*m_w+j].m_color.m_rgb[1]=color.m_rgb[1]*mag;
		m_pixel[i*m_w+j].m_color.m_rgb[2]=color.m_rgb[2]*mag;
	/*	if (m_pixel[i*m_w+j].m_color.m_rgb[0] > 0 ||
			m_pixel[i*m_w+j].m_color.m_rgb[1] > 0|| 
			m_pixel[i*m_w+j].m_color.m_rgb[2] > 0) {

			printf("non-zero");
		}*/
	}

private:


	//////////////////////////////////////////////
	// over_sample
	/////////////////////////////////////////////
	bool OverSampler(Vector3f* m_sample, int over_sample_ratio){
		float u_step=1.0f/(m_w*over_sample_ratio);
		float v_step=1.0f/(m_h*over_sample_ratio);

		float u=u_step/2.0f, v=v_step/2.0f;
		// cout<<LL<<UL<<LR<<UR<<endl;
		for (int i=0; i<m_h*over_sample_ratio; i++){
			for(int j=0; j< m_w*over_sample_ratio; j++){
				m_sample[i*m_w*over_sample_ratio+j]=u*(v*UR+(1-v)*LR)+(1-u)*(v*UL+(1-v)*LL);
				// m_pixel[i*m_w+j].m_color.m_rgb[0]=0.0;
				// m_pixel[i*m_w+j].m_color.m_rgb[1]=0.0;
				// m_pixel[i*m_w+j].m_color.m_rgb[2]=0.0;
				// cout<<m_sample[i*m_w*over_sample_ratio+j]<<endl;
				u+=u_step;

			}
			// cout<<"u: "<<u << ", "<<v<<endl;
			u=u_step/2.0f;
			v+=v_step;
		}
		return true;
	}


	//////////////////////////////////////////////
	// jitter_over_sample
	/////////////////////////////////////////////
	float Rand(float step){
		return ((((float) rand() / RAND_MAX)-.5f)*step); //bound
	}
	bool JitterSampler(Vector3f* m_sample, int over_sample_ratio){

		float u_step=1.0f/m_w/over_sample_ratio;
		float v_step=1.0f/m_h/over_sample_ratio;

		float u=u_step/2.0f, v=v_step/2.0f;
		float u_r=u+Rand(u_step), v_r=v+Rand(v_step);
		// cout<<LL<<UL<<LR<<UR<<endl;
		for (int i=0; i<m_h*over_sample_ratio; i++){
			for(int j=0; j< m_w*over_sample_ratio; j++){
				m_sample[i*m_w*over_sample_ratio+j]=u_r*(v_r*UR+(1-v_r)*LR)+(1-u_r)*(v_r*UL+(1-v_r)*LL);
				// m_pixel[i*m_w+j].m_color.m_rgb[0]=0.0;
				// m_pixel[i*m_w+j].m_color.m_rgb[1]=0.0;
				// m_pixel[i*m_w+j].m_color.m_rgb[2]=0.0;
				// cout<<m_sample[i*m_w*over_sample_ratio+j]<<endl;
				u+=u_step;
				u_r=u+Rand(u_step);
				v_r=v+Rand(v_step);
				// cout<<u_step<<", "<<v_step<<endl;

			}
			// cout<<"u: "<<u << ", "<<v<<endl;
			u=u_step/2.0f;
			v+=v_step;
		}
		return true;
	}



public:

	bool Sample(int depth, int over_sample_ratio, SamplerType s){
		m_sample=new Vector3f[m_w*m_h*over_sample_ratio*over_sample_ratio];
		if (s==OverS){
			OverSampler(m_sample, over_sample_ratio);
		}
		else if (s==JitterS){
			JitterSampler(m_sample, over_sample_ratio);
		}
		else {
			cerr<<"wrong sampler type, either: OverS or JitterS"<<endl;
			return false;
		}
		
		//float t_min, t_max;
		//CRayTracer ray_tracer;
		float num_samples = (float) over_sample_ratio * over_sample_ratio; 
		//#pragma omp parallel for
		for(int i=0; i< m_h; i++) {
			for(int j=0; j<m_w; j++){
			/*	if (i != 150 || j != 101)
					continue; */
				CColor temp(0.0f, 0.0f, 0.0f);
				Vector3f pos, dir;
				for(int k=i; k<i+over_sample_ratio;k++) {
					for(int g=j;g<j+over_sample_ratio;g++){
						pos=m_sample[k*m_h*over_sample_ratio+g];
						dir=pos-m_eye;
						
						//CRay ray(pos, dir, dir.norm(), INF);
						CRay ray(pos, dir, EPS, INF);
						CColor ray_color;
						CBRDF brdf; 
						bool isObj = m_rayTracer->trace(ray, depth, &ray_color, brdf);  //check ???? if ray_tracer behave right???
						temp = temp.Add(ray_color);		///completed cumulated oversampling
						if (isObj) 
							temp = temp.Add(brdf.ka);
					}
					
					//temp.Print(); 
					
					
				}

				//temp.Print();
				temp.m_rgb[0] /=num_samples;
				temp.m_rgb[1] /= num_samples;
				temp.m_rgb[2] /= num_samples;
				/*if (temp.m_rgb[0] > 0.0 ||
					temp.m_rgb[1] > 0.0 ||
					temp.m_rgb[2] > 0.0) {
					temp.Print();
				}*/

				ColorPixel(i,j,temp);

			/*	m_pixel[i*m_w+j].m_color =  m_pixel[i*m_w+j].m_color.Add(temp);
			
				m_pixel[i*m_w+j].m_color.m_rgb[0] /= num_samples;
				m_pixel[i*m_w+j].m_color.m_rgb[1] /= num_samples;
				m_pixel[i*m_w+j].m_color.m_rgb[2] /= num_samples;*/
					// m_pixel[i*m_w+j].m_color.Print();
			}
		}
			
		return true;
	}

	bool Film(){
		CImg<unsigned char> img(m_w, m_h, 1, 3);
		img.fill(0);
	/*	m_pixel[200].m_color.m_rgb[0] = 1;
		m_pixel[200].m_color.m_rgb[1] = 1;
	    m_pixel[200].m_color.m_rgb[2] = 1;*/ 
		cimg_forXYC(img, x, y, c) {img(x,y,c)=(unsigned char) m_pixel[x*m_w+y].m_color.m_rgb[c];}
		//img.normalize(0,255);
		img.save(g_fname.c_str());
		img.display("RayTracer");
		return true;
	}



	void DestroyCamera(){
		delete[] m_pixel;
		delete[] m_sample;
	}


public:
	int m_w, m_h; // width and height  in pixels
	float m_width, m_height; //width and height in world units (cm)
	// float m_u_step, m_v_step, m_u, m_v;  //parameters base one (0,1) for width (u), height (v)
	Vector3f m_center;
	Vector3f m_normal;
	Vector3f LL, UL, LR, UR;
	Vector3f m_eye;
	// int m_over_sample_ratio;
	Pixel* m_pixel;
	Vector3f* m_sample;
	CRayTracer* m_rayTracer; // = new CRayTracer(); 
};

<<<<<<< HEAD
CPrimitive* InitScene() { 
	vector<CPrimitive*> primList; 
	CBRDF brdf; 

	brdf.ka = CColor(0.0f, 0.0f, 0.0f);   // ambient 
	brdf.kd = CColor(1.0f, 0.0f, 0.0f);   // diffuse
	brdf.ks = CColor(0.0f, 0.0f, 0.0f);   // specular
	//brdf.ks = CColor(0.2f, 0.2f, 0.2f);   // specular
	brdf.kr = CColor(1.0f, 1.0f, 1.0f);   // reflection
	brdf.p  = 16.0f;         
	//CMaterial* mat = ;

	CMatrix m;
	// m.scale(Vector3f(1,2,1));
	// m.translate(Vector3f(0,0,-10));

	// CTransformation T(m); 
	CGeometricPrimitive* prim1 = new CGeometricPrimitive(); 
	prim1->m_objToWorld = CTransformation(m); 
	
	CMatrix M;
	// M.translate(Vector3f(0,0,10));

	prim1->m_worldToObj = CTransformation(M); 
	prim1->m_mat = new CMaterial(brdf); 
	prim1->m_shape = new CSphere(V3f(0.0f, 0.0f, 0.0f), 5.0f);// new CTriangle(V3f(0.0f, 0.0f, -5.0f), V3f(5.0, 0.0, -5.0f), V3f(0.0, 5.0f, -5.0f));//new CSphere(V3f(-6.0f, 0.0f, -5.0f), 5.0f);
	primList.push_back(prim1);
	// CGeometricPrimitive* prim2 = new CGeometricPrimitive(); 
	// prim2->m_objToWorld = T; prim2->m_worldToObj = T; 
	// brdf.kd =  CColor(0.0f, 1.0f, 0.0f);   // diffuse
	// prim2->m_mat = new CMaterial(brdf); 
	// prim2->m_shape = new CSphere(V3f(5.0f, -0.0f, -5.0f), 5.0f);// new CTriangle(V3f(0.0f, 0.0f, -10.0f), V3f(15.0, 0.0, -10.0f), V3f(0.0, 15.0f, -10.0f));
	// primList.push_back(prim2);

	CAggregatePrimitive* scene = new CAggregatePrimitive(primList);
	//scene->m_objToWorld = CTransformation();
	//scene->m_worldToObj = CTransformation();
	              // specular 
	//scene->m_mat = new CMaterial(brdf);
	//scene->m_shape = new CSphere(V3f(0.0f, 0.0f, -5.0f), 5.0f);

	//scene->m_shape = new CTriangle(V3f(0.0f, 0.0f, -5.0f), V3f(15.0, 0.0, -0.0f), V3f(0.0, 15.0f, -5.0f));
	return (CPrimitive*)scene; 	
}


void ParseObject() {

}

vector<CLight*> InitLights() {
	vector<CLight*> lights; 
	CLight* light1 = new CLight(V3f(0.0f, 0.0f, 50.0f), CColor(1.0f, 1.0f, 1.0f), CLight::Point);
	lights.push_back(light1);
	return lights; 
=======
//CPrimitive* InitScene() { 
//	vector<CPrimitive*> primList; 
//	CBRDF brdf; 
//
//	brdf.ka = CColor(1.0f, 1.0f, 1.0f);   // ambient 
//	brdf.kd = CColor(0.0f, 0.0f, 0.0f);   // diffuse
//	brdf.ks = CColor(0.0f, 0.0f, 0.0f);   // specular
//	//brdf.ks = CColor(0.2f, 0.2f, 0.2f);   // specular
//	brdf.kr = CColor(0.0f, 0.0f, 0.0f);   // reflection
//	brdf.p  = 16.0f;         
//	//CMaterial* mat = ;
//	CTransformation T; 
//	CGeometricPrimitive* prim1 = new CGeometricPrimitive(); 
//	prim1->m_objToWorld = T; 
//	prim1->m_worldToObj = T; 
//	prim1->m_mat = new CMaterial(brdf); 
//	prim1->m_shape = new CSphere(V3f(0.0f, -0.0f, -5.0f), 5.0f);// new CTriangle(V3f(0.0f, 0.0f, -5.0f), V3f(5.0, 0.0, -5.0f), V3f(0.0, 5.0f, -5.0f));//new CSphere(V3f(-6.0f, 0.0f, -5.0f), 5.0f);
//	primList.push_back(prim1);
//	//CGeometricPrimitive* prim2 = new CGeometricPrimitive(); 
//	//prim2->m_objToWorld = T; prim2->m_worldToObj = T; 
//	//brdf.kd =  CColor(0.0f, 1.0f, 0.0f);   // diffuse
//	//prim2->m_mat = new CMaterial(brdf); 
//	//prim2->m_shape = new CSphere(V3f(5.0f, -0.0f, -5.0f), 5.0f);// new CTriangle(V3f(0.0f, 0.0f, -10.0f), V3f(15.0, 0.0, -10.0f), V3f(0.0, 15.0f, -10.0f));
//	//primList.push_back(prim2);
//
//	CAggregatePrimitive* scene = new CAggregatePrimitive(primList);
//	//scene->m_objToWorld = CTransformation();
//	//scene->m_worldToObj = CTransformation();
//	              // specular 
//	//scene->m_mat = new CMaterial(brdf);
//	//scene->m_shape = new CSphere(V3f(0.0f, 0.0f, -5.0f), 5.0f);
//
//	//scene->m_shape = new CTriangle(V3f(0.0f, 0.0f, -5.0f), V3f(15.0, 0.0, -0.0f), V3f(0.0, 15.0f, -5.0f));
//	return (CPrimitive*)scene; 	
//}


//vector<CLight*> InitLights() {
//	vector<CLight*> lights; 
//	CLight* light1 = new CLight(V3f(0.0f, 100.0f, 0.0f), CColor(1.0f, 1.0f, 1.0f), CLight::Point);
//	lights.push_back(light1);
//	return lights; 
//}




void loadScene(string file) {
  //store variables and set stuff at the end
  g_lights.clear(); 
  vector<CPrimitive*> prims; 
  vector<V3f> vertices; 
  vertices.clear(); 
  CBRDF brdf; 
  CTransformation M; 
  g_fname = "output.bmp";

  std::ifstream inpfile(file.c_str());
  if(!inpfile.is_open()) {
    std::cout << "Unable to open file" << std::endl;
  } else {
    std::string line;
    //MatrixStack mst;

    while(inpfile.good()) {
      std::vector<std::string> splitline;
      std::string buf;
	   
      std::getline(inpfile,line);
      std::stringstream ss(line);

      while (ss >> buf) {
        splitline.push_back(buf);
      }
      //Ignore blank lines
      if(splitline.size() == 0) {
        continue;
      }

      //Ignore comments
      if(splitline[0][0] == '#') {
        continue;
      }

      //Valid commands:
      //size width height
      //  must be first command of file, controls image size
      else if(!splitline[0].compare("size")) {
        m_width = atoi(splitline[1].c_str());
        m_height = atoi(splitline[2].c_str());
      }
      //maxdepth depth
      //  max # of bounces for ray (default 5)
      else if(!splitline[0].compare("maxdepth")) {
        // maxdepth: atoi(splitline[1].c_str())
      }
      //output filename
      //  output file to write image to 
      else if(!splitline[0].compare("output")) {
        g_fname = splitline[1];
      }

      //camera lookfromx lookfromy lookfromz lookatx lookaty lookatz upx upy upz fov
      //  specifies the camera in the standard way, as in homework 2.
      else if(!splitline[0].compare("camera")) {
        // lookfrom:
        //    atof(splitline[1].c_str())
        //    atof(splitline[2].c_str())
        //    atof(splitline[3].c_str())
        // lookat:
        //    atof(splitline[4].c_str())
        //    atof(splitline[5].c_str())
        //    atof(splitline[6].c_str())
        // up:
        //    atof(splitline[7].c_str())
        //    atof(splitline[8].c_str())
        //    atof(splitline[9].c_str())
        // fov: atof(splitline[10].c_str());
      }

      //sphere x y z radius
      //  Deﬁnes a sphere with a given position and radius.
      else if(!splitline[0].compare("sphere")) {
		 CGeometricPrimitive* sphere = new CGeometricPrimitive(); 
		 sphere->m_worldToObj = M; 
		 sphere->m_objToWorld = M; // BUG
         float x = (float)atof(splitline[1].c_str());
         float y = (float)atof(splitline[2].c_str());
         float z = (float)atof(splitline[3].c_str());
         float r = (float)atof(splitline[4].c_str());
		 sphere->m_mat = new CMaterial(brdf);
		 sphere->m_shape = new CSphere(V3f(x, y, z-5), r);
		 prims.push_back(sphere);
        // Create new sphere:
        //   Store 4 numbers
        //   Store current property values
        //   Store current top of matrix stack
      }
      //maxverts number
      //  Deﬁnes a maximum number of vertices for later triangle speciﬁcations. 
      //  It must be set before vertices are deﬁned.
      else if(!splitline[0].compare("maxverts")) {
        // Care if you want
        // Here, either declare array size
        // Or you can just use a STL vector, in which case you can ignore this
      }
      //maxvertnorms number
      //  Deﬁnes a maximum number of vertices with normals for later speciﬁcations.
      //  It must be set before vertices with normals are deﬁned.
      else if(!splitline[0].compare("maxvertnorms")) {
        // Care if you want
      }
      //vertex x y z
      //  Deﬁnes a vertex at the given location.
      //  The vertex is put into a pile, starting to be numbered at 0.
      else if(!splitline[0].compare("vertex")) {

        float x = (float)atof(splitline[1].c_str());
        float y = (float)atof(splitline[2].c_str());
        float z = (float)atof(splitline[3].c_str());
		vertices.push_back(V3f(x, y, z-3));

        // Create a new vertex with these 3 values, store in some array
      }
      //vertexnormal x y z nx ny nz
      //  Similar to the above, but deﬁne a surface normal with each vertex.
      //  The vertex and vertexnormal set of vertices are completely independent
      //  (as are maxverts and maxvertnorms).
      else if(!splitline[0].compare("vertexnormal")) {
        // x: atof(splitline[1].c_str()),
        // y: atof(splitline[2].c_str()),
        // z: atof(splitline[3].c_str()));
        // nx: atof(splitline[4].c_str()),
        // ny: atof(splitline[5].c_str()),
        // nz: atof(splitline[6].c_str()));
        // Create a new vertex+normal with these 6 values, store in some array
      }
      //tri v1 v2 v3
      //  Create a triangle out of the vertices involved (which have previously been speciﬁed with
      //  the vertex command). The vertices are assumed to be speciﬁed in counter-clockwise order. Your code
      //  should internally compute a face normal for this triangle.
      else if(!splitline[0].compare("tri")) {
        int i = (int)atof(splitline[1].c_str());
        int j = (int)atof(splitline[2].c_str());
        int k = (int)atof(splitline[3].c_str());
		CGeometricPrimitive* tri = new CGeometricPrimitive(); 
		tri->m_worldToObj = M; 
		tri->m_objToWorld = M; // BUG

		tri->m_mat = new CMaterial(brdf);
		tri->m_shape = new CTriangle(vertices[i], vertices[j], vertices[k]);
		prims.push_back(tri);
        // Create new triangle:
        //   Store pointer to array of vertices
        //   Store 3 integers to index into array
        //   Store current property values
        //   Store current top of matrix stack
      }
      //trinormal v1 v2 v3
      //  Same as above but for vertices speciﬁed with normals.
      //  In this case, each vertex has an associated normal, 
      //  and when doing shading, you should interpolate the normals 
      //  for intermediate points on the triangle.
      else if(!splitline[0].compare("trinormal")) {
        // v1: atof(splitline[1].c_str())
        // v2: atof(splitline[2].c_str())
        // v3: atof(splitline[3].c_str())
        // Create new triangle:
        //   Store pointer to array of vertices (Different array than above)
        //   Store 3 integers to index into array
        //   Store current property values
        //   Store current top of matrix stack
      }

      //translate x y z
      //  A translation 3-vector
      else if(!splitline[0].compare("translate")) { //TO-DO
        float x = atof(splitline[1].c_str()); 
        float y = atof(splitline[2].c_str());
        float z = atof(splitline[3].c_str());
		
        // Update top of matrix stack
      }
      //rotate x y z angle
      //  Rotate by angle (in degrees) about the given axis as in OpenGL.
      else if(!splitline[0].compare("rotate")) { //TO-DO 
        float x = (float)atof(splitline[1].c_str());
        float y = (float)atof(splitline[2].c_str());
        float z = (float)atof(splitline[3].c_str());
        float angle = (float)atof(splitline[4].c_str());
        // Update top of matrix stack
      }
      //scale x y z
      //  Scale by the corresponding amount in each axis (a non-uniform scaling).
      else if(!splitline[0].compare("scale")) {
        float x = (float)atof(splitline[1].c_str());
        float y = (float)atof(splitline[2].c_str());
        float z = (float)atof(splitline[3].c_str());
        // Update top of matrix stack
      }
      //pushTransform
      //  Push the current modeling transform on the stack as in OpenGL. 
      //  You might want to do pushTransform immediately after setting 
      //   the camera to preserve the “identity” transformation.
      else if(!splitline[0].compare("pushTransform")) {
        //mst.push();
      }
      //popTransform
      //  Pop the current transform from the stack as in OpenGL. 
      //  The sequence of popTransform and pushTransform can be used if 
      //  desired before every primitive to reset the transformation 
      //  (assuming the initial camera transformation is on the stack as 
      //  discussed above).
      else if(!splitline[0].compare("popTransform")) {
        //mst.pop();
      }

      //directional x y z r g b
      //  The direction to the light source, and the color, as in OpenGL.
      else if(!splitline[0].compare("directional")) {
        float x = (float)atof(splitline[1].c_str()); 
        float y = (float)atof(splitline[2].c_str());
        float z = (float)atof(splitline[3].c_str());
        float r = (float)atof(splitline[4].c_str());
        float g = (float)atof(splitline[5].c_str());
        float b = (float)atof(splitline[6].c_str());
		CLight* light = new CLight(V3f(x, y, z), CColor(r, g, b), CLight::Directional);
		g_lights.push_back(light);
        // add light to scene...
      }
      //point x y z r g b
      //  The location of a point source and the color, as in OpenGL.
      else if(!splitline[0].compare("point")) {
		  float x = (float)atof(splitline[1].c_str()); 
		  float y = (float)atof(splitline[2].c_str());
		  float z = (float)atof(splitline[3].c_str());
		  float r = (float)atof(splitline[4].c_str());
		  float g = (float)atof(splitline[5].c_str());
		  float b = (float)atof(splitline[6].c_str());
		  CLight* light = new CLight(V3f(x, y, z), CColor(r, g, b), CLight::Point);
		  g_lights.push_back(light);
        // add light to scene...
      }
      //attenuation const linear quadratic
      //  Sets the constant, linear and quadratic attenuations 
      //  (default 1,0,0) as in OpenGL.
      else if(!splitline[0].compare("attenuation")) {
        // const: atof(splitline[1].c_str())
        // linear: atof(splitline[2].c_str())
        // quadratic: atof(splitline[3].c_str())
      }
      //ambient r g b
      //  The global ambient color to be added for each object 
      //  (default is .2,.2,.2)
      else if(!splitline[0].compare("ambient")) {
        float r = (float)atof(splitline[1].c_str());
        float g = (float)atof(splitline[2].c_str());
        float b = (float)atof(splitline[3].c_str());
		brdf.ka =CColor(r,g,b);
      }

      //diﬀuse r g b
      //  speciﬁes the diﬀuse color of the surface.
      else if(!splitline[0].compare("diffuse")) {
		  float r = (float)atof(splitline[1].c_str());
		  float g = (float)atof(splitline[2].c_str());
		  float b = (float)atof(splitline[3].c_str());
		  brdf.kd = CColor(r,g,b);
        // Update current properties
      }
      //specular r g b 
      //  speciﬁes the specular color of the surface.
      else if(!splitline[0].compare("specular")) {
		  float r = (float)atof(splitline[1].c_str());
		  float g = (float)atof(splitline[2].c_str());
		  float b = (float)atof(splitline[3].c_str());
		  brdf.ks = CColor(r,g,b);
        // Update current properties
      }
      //shininess s
      //  speciﬁes the shininess of the surface.
      else if(!splitline[0].compare("shininess")) {
        float p = (float)atof(splitline[1].c_str());
		brdf.p = p; 
        // Update current properties
      }
      //emission r g b
      //  gives the emissive color of the surface.
      else if(!splitline[0].compare("emission")) {
        // r: atof(splitline[1].c_str())
        // g: atof(splitline[2].c_str())
        // b: atof(splitline[3].c_str())
        // Update current properties
      } else {
        std::cerr << "Unknown command: " << splitline[0] << std::endl;
      }
    }

    inpfile.close();
  }

  vertices.clear(); 
  g_scene = new CAggregatePrimitive(prims); 
>>>>>>> FETCH_HEAD
}

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("invalid input\n");
		exit(-1);
	}

	loadScene(string(argv[1]));


<<<<<<< HEAD
	Vector3f eye(0,0,5);
	int w = 100; 
	int h = w;
=======
	Vector3f eye(0,0,10);
	//int w = 1024; 
	//int h = w;
>>>>>>> FETCH_HEAD
	Vector3f LL(-10, -10, 0), UL(-10,10, 0),
		LR(10, -10, 0), UR(10, 10, 0);
	int over_sample_ratio = 1;
	CColor pixel(1,0,0);

	CCamera camera(eye, m_width, m_height, LL, UL, LR, UR);
	CRayTracer* rayTracer = new CRayTracer(); 
<<<<<<< HEAD
	vector<CLight*> lights = InitLights(); 

	CPrimitive* scene = InitScene(); 

	rayTracer->Setup(scene, lights);

=======
	//vector<CLight*> lights = InitLights(); 
	//CPrimitive* scene = InitScene(); 
	 
	
	rayTracer->Setup(g_scene, g_lights);
>>>>>>> FETCH_HEAD
	camera.SetupRayTracer(rayTracer);

	// omp_set_num_threads(12);

	
	camera.Sample(1, over_sample_ratio, CCamera::OverS);
	//DELETE_OBJECT(timer);


	camera.Film();

	FOR_u (i, g_lights.size())
		DELETE_OBJECT(g_lights[i]);

	DELETE_OBJECT(g_scene);
	DELETE_OBJECT(rayTracer);
	return 0;

}

