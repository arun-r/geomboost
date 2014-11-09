#include <geomboost.h>
#include <vector>
#include <iostream>
#include <time.h>

#include <map>
#include <list>
#include <Windows.h>

#undef max


static const float PI = 3.14159265358979323846f;

struct testpoint2{
	float xl,yl,zl;
	float x()const{
		return xl;
	}
	float y()const{
		return yl;
	}
	float z()const{
		return zl;
	}
};


struct testpoint{
	float xl,yl,zl;
	testpoint():xl(0),yl(0),zl(0){}
	testpoint(float x1, float y1, float z1):
		xl(x1), yl(y1), zl(z1){}

	float x()const{
		return xl;
	}
	float y()const{
		return yl;
	}
	float z()const{
		return zl;
	}
};

std::ostream& operator<<(std::ostream& ostr, const testpoint& tp) {
	ostr << tp.x() << ":" << tp.y() << ":" << tp.z();
	return ostr;
}

using namespace std;
using namespace geomboost;

static void test_on_sphere(float avgradius, float thickness, int numtests){
	std::vector<testpoint> vp;
	srand ((unsigned int)5);

	const float numpoints = 1000000.0f;
	
	float jjincr = PI * sqrt(1.0f/numpoints);
	float iiincr = 2*PI * sqrt(1.0f/numpoints);

	

	for(float ii = 0; ii < 2*PI; ii += iiincr){
		for (float jj = 0; jj < PI; jj += jjincr){
			float radius = ((float)rand()/RAND_MAX - 0.5f)*thickness + avgradius;
			float cosii = cos(ii);
			float sinii = sin(ii);
			float cosjj = cos(jj);
			float sinjj = sin(jj);
			vp.push_back(testpoint(radius*cosii*sinjj, radius*sinii*sinjj, radius*cosjj));
		}
	}
	closest_point_finder<testpoint> cpf(vp.begin(), vp.end()); 
	
	float totaltime1 = 0;
	float totaltime2 = 0;
	for(int ii = 0; ii < numtests; ++ii)
	{
		testpoint tp(((float)rand()/RAND_MAX - 0.5f)*thickness*1.5f + avgradius, 
					((float)rand()/RAND_MAX - 0.5f)*thickness*1.5f + avgradius, 
					((float)rand()/RAND_MAX - 0.5f)*thickness*1.5f + avgradius);
		//testpoint tp((float)(8+rand()%4), (float)(8+rand()%4), (float)(8+rand()%4));
		DWORD t1, t2, t3;
		t1 = GetTickCount();
		const testpoint* cp1 = cpf.find_closest(tp);
		t2 = GetTickCount();
		const testpoint* cp2 = cpf.find_closest_linear(tp);
		t3 = GetTickCount();
		totaltime1 += (t2 -t1);
		totaltime2 += (t3 -t2);
		float d1 = geomboost::distance(*cp1, tp);
		float d2 = geomboost::distance(*cp2, tp);
		if( fabs(d1 - d2) > 0)
		{	
			cout << " AT " << ii << " dist1 = " << geomboost::distance(*cp1, tp) << ":" << " TIME = " << t2-t1;
			cout << " dist2 = " << geomboost::distance(*cp2, tp) << ":" << " TIME = " << t3-t2 << std::endl;			;
			cout << " TEST POINT = " << tp << " D1 = " << *cp1 << " D2 = " << *cp2 << endl;			
		}
		if (ii % 10 == 0){
			cout << " AVG TIME1 = " << totaltime1/(ii+1) << "TIME 2 " << totaltime2/(ii+1) << endl;
		}
	}

	cout << " AVG TIME1 = " << totaltime1/numtests << "TIME 2 " << totaltime2/numtests << endl;
}



void tester_pod_xyz(){

	
	std::map<int, testpoint> mp;
	std::map<testpoint, int> mp2;
	std::list<const testpoint*> lpp;
	std::list<testpoint*> lpp2;
	std::list<const testpoint> lpp3;

	
	/*

	closest_point_finder<testpoint> cpf(vp.begin(), vp.end());

	
	closest_point_finder<testpoint> cpf2(mp.begin(), mp.end());
	closest_point_finder<testpoint> cpf3(mp2.begin(), mp2.end());
	closest_point_finder<testpoint> cpf4(lpp.begin(), lpp.end());
	closest_point_finder<testpoint> cpf5(lpp2.begin(), lpp2.end());
	closest_point_finder<testpoint> cpf6(lpp3.begin(), lpp3.end());
	closest_point_finder<testpoint> cpf7(pps, pps+5);

	const testpoint* cp = cpf.find_closest(pps[3]);
	const testpoint* cp1 = cpf2.find_closest(pps[3]);
	const testpoint* cp3 = cpf3.find_closest(pps[3]);
	const testpoint* cp4 = cpf4.find_closest(pps[3]);
	const testpoint* cp5 = cpf5.find_closest(pps[3]);
	const testpoint* cp6 = cpf6.find_closest(pps[3]);
	const testpoint* cp7 = cpf7.find_closest(pps[3]);


	if (cp)	std::cout << *cp << " dist " << distance(*cp, pps[3]) << std::endl;
	if (cp1)	std::cout << *cp1;
	if (cp3)	std::cout << *cp3;
	if (cp4)	std::cout << *cp4;
	if (cp5)	std::cout << *cp5;
	if (cp6)	std::cout << *cp6;
	if (cp7)	std::cout << *cp7;

	

	boundingbox<testpoint> bb(lpp.begin(), lpp.end());
	boundingbox<testpoint> bb2(lpp2.begin(), lpp2.end());
	boundingbox<testpoint> bb3(lpp3.begin(), lpp3.end());

	testpoint t1, t2;
	testpoint t = t1 + t2;*/

	


} 

int main(){
	test_on_sphere(5, 2, 100);
	return 0;

	 
}
