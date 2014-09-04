#include <geomboost.h>
#include <vector>
#include <iostream>

#include <map>
#include <list>


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
	ostr << tp.x() << ":" << tp.y() << ":" << tp.z() << std::endl;
	return ostr;
}


using namespace geomboost;



void tester_pod_xyz(){

	std::vector<testpoint> vp;
	std::map<int, testpoint> mp;
	std::map<testpoint, int> mp2;
	std::list<const testpoint*> lpp;
	std::list<testpoint*> lpp2;
	std::list<const testpoint> lpp3;

	float a = 5,b = 10,c = 5;

	const unsigned int ii = 10;
	

	unsigned kk = std::max((unsigned int)((a-b)/c), ii);

	testpoint pps[10];
	for(float ii = -5; ii < 5; ++ii){
		vp.push_back(testpoint(ii,ii+1,ii-1));
	}

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


	if (cp)	std::cout << *cp;
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
	testpoint t = t1 + t2;

	

}

int main(){
	tester_pod_xyz();
	return 0;
}
