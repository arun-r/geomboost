
#ifndef SPACE_HASH_H
#define SPACE_HASH_H


#include <map>
#include <vector>
#include <list>
#include <type_traits>
#include <algorithm>

#include <boost/utility/declval.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits.hpp>


namespace geomboost{


#define HAS_MEMBER_FUNC(structname, funcname)						\
	template <typename T, typename FuncSign>						\
	struct structname{												\
		typedef char yes[1];										\
		typedef char no [2];										\
		template <typename U, U> struct type_check;					\
		template <typename _1> static yes &chk(type_check<FuncSign, &_1::funcname > *);	\
		template <typename   > static no  &chk(...);							\
		static bool const value = (sizeof(chk<T>(0)) == sizeof(yes));			\
	}

HAS_MEMBER_FUNC(has_x, x);
HAS_MEMBER_FUNC(has_y, y);
HAS_MEMBER_FUNC(has_z, z);



template <typename T>
struct t_traits{
	typedef BOOST_TYPEOF_TPL(boost::declval<T const&>().x()) real_t;

	static real_t getx(const T& val){
		return val.x();
	}

	static real_t gety(const T& val){
		return val.y();
	}

	typedef real_t (T::*ZFunc)()const;

	template <typename T1>
	static real_t getz(const T1& val, typename std::enable_if<has_z<T1, ZFunc>::value>::type* = 0){
		return val.z();
	}

	template <typename T1>
	static real_t getz(const T1& val, typename std::enable_if<!has_z<T1, ZFunc>::value>::type* = 0){
		return 0;
	}

	template <typename T1 >
	static T1 makeT(real_t x, real_t y, real_t z, typename std::enable_if<boost::is_pod<T1>::value>::type* = 0){
		T1 tobj = {x,y,z};
		return tobj;
	}

	template <typename T1 >
	static T1 makeT(real_t x, real_t y, real_t z, typename std::enable_if<!boost::is_pod<T1>::value>::type* = 0){
		return T1 (x,y,z);
	}

	template <typename T1 >
	static T1 makeT(real_t x, real_t y, typename std::enable_if<boost::is_pod<T1>::value>::type* = 0){
		T1 tobj = {x,y};
		return tobj;
	}

	template <typename T1 >
	static T1 makeT(real_t x, real_t y, typename std::enable_if<!boost::is_pod<T1>::value>::type* = 0){
		return T1(x,y);
	}


	template <typename T1 >
	static T1 sum(const T1& lhs, const T1& rhs, typename std::enable_if<has_z<T1, ZFunc>::value>::type* = 0){	
		return makeT<T1>(getx(lhs) + getx(rhs), gety(lhs) + gety(rhs), getz(lhs) + getz(rhs));
	}

	template <typename T1 >
	static T1 sum(const T1& lhs, const T1& rhs, typename std::enable_if<!has_z<T1, ZFunc>::value>::type* = 0){	
		return makeT<T1>(getx(lhs) + getx(rhs), gety(lhs) + gety(rhs));
	}

	template <typename T1 >
	static T1 diff(const T1& lhs, const T1& rhs, typename std::enable_if<has_z<T1, ZFunc>::value>::type* = 0){	
		return makeT<T1>(getx(lhs) - getx(rhs),	gety(lhs) - gety(rhs), getz(lhs) - getz(rhs));
	}

	template <typename T1 >
	static T1 diff(const T1& lhs, const T1& rhs, typename std::enable_if<!has_z<T1, ZFunc>::value>::type* = 0){	
		return makeT<T1>(getx(lhs) - getx(rhs), gety(lhs) - gety(rhs));
	}

	template <typename T1>
	static T1 cp(const T1& lhs, const T1& rhs, typename std::enable_if<has_z<T1, ZFunc>::value>::type* = 0){
		return makeT<T1>(gety(lhs)*getz(rhs) - getz(lhs)*gety(rhs), 
					-(getx(lhs)*getz(rhs) - getx(rhs)*getz(lhs)),
					getx(lhs)*gety(rhs) - getx(rhs)*gety(lhs));
	}

	template <typename T1>
	static T1 cp(const T1& lhs, const T1& rhs, typename std::enable_if<!has_z<T1, ZFunc>::value>::type* = 0){
		throw "bad call";
	}
};



template <typename T>
T operator+(const T& lhs, const T& rhs){
	return t_traits<T>::sum(lhs, rhs);
}

template <typename T>
T operator-(const T& lhs, const T& rhs){
	return t_traits<T>::diff(lhs, rhs, true);
}

template <typename T>
typename t_traits<T>::real_t squared_distance(const T& lhs, const T& rhs){
	t_traits<T>::real_t xdiff = t_traits<T>::getx(lhs) - t_traits<T>::getx(rhs);
	t_traits<T>::real_t ydiff = t_traits<T>::gety(lhs) - t_traits<T>::gety(rhs);
	t_traits<T>::real_t zdiff = t_traits<T>::getz(lhs) - t_traits<T>::getz(rhs);
	return xdiff*xdiff + ydiff*ydiff + zdiff*zdiff;
}

template <typename T>
typename t_traits<T>::real_t distance(const T& lhs, const T& rhs){
	return sqrt(squared_distance(lhs, rhs));
}

template <typename T>
typename t_traits<T>::real_t magnitude(const T& val){
	t_traits<T>::real_t xval = t_traits::getx(val);
	t_traits<T>::real_t yval = t_traits::gety(val);
	t_traits<T>::real_t zval = t_traits::getz(val);
	return xval*xval + yval*yval + zval*zval;
}

template <typename T>
typename t_traits<T>::real_t dotproduct(const T& lhs, const T& rhs){
	return t_traits::getx(lhs)*t_traits::getx(rhs) +
			t_traits::gety(lhs)*t_traits::gety(rhs) +
			t_traits::getz(lhs)*t_traits::getz(rhs);
}

template <typename T>
T crossproduct(const T& lhs, const T& rhs){
	return t_traits::cp(lhs, rhs);
}



template <typename T>
struct boundingbox : public t_traits<T>{
private:
	void init(const std::vector<const T*>& points){

		if(points.size() == 0){
			xmin = xmax = ymin = ymax = zmin = zmax = 0;
			return;
		}

		xmin = xmax = getx(*points.front());
		ymin = ymax = gety(*points.front());
		zmin = zmax = getz(*points.front());

		for(std::vector<const T*>::const_iterator iter = points.begin(); 
				iter != points.end(); ++iter){
		
			xmin = std::min(xmin, getx(**iter));
			ymin = std::min(ymin, gety(**iter));
			zmin = std::min(zmin, getz(**iter));

			xmax = std::max(xmax, getx(**iter));
			ymax = std::max(ymax, gety(**iter));
			zmax = std::max(zmax, getz(**iter));
		}
	}

public:
	real_t xmin, xmax;
	real_t ymin, ymax;
	real_t zmin, zmax;

	template <typename IterT>
	boundingbox(IterT begin, IterT end,
		typename std::enable_if< std::is_same<typename std::iterator_traits<IterT>::value_type, T >::value, void >::type* dummy = 0){
		std::vector<const T*> tvec;
		for(IterT iter = begin; iter != end; ++iter){
			tvec.push_back(&(*iter));
		}
		init(tvec);
	}

	template <typename IterT>
	boundingbox(IterT begin, IterT end,
		typename std::enable_if< std::is_same<typename std::iterator_traits<IterT>::value_type, T* >::value>::type* dummy = 0){		
		std::vector<const T*> tvec;
		for(IterT iter = begin; iter != end; ++iter){
			tvec.push_back((*iter));
		}
		init(tvec);
	}

	template <typename IterT>
	boundingbox(IterT begin, IterT end,
		typename std::enable_if< std::is_same<typename std::iterator_traits<IterT>::value_type, const T* >::value>::type* dummy = 0){
		
		std::vector<const T*> tvec;
		for(IterT iter = begin; iter != end; ++iter){
			tvec.push_back((*iter));
		}
		init(tvec);
	}

	boundingbox():xmin(0), xmax(0), ymin(0), ymax(0), zmin(0), zmax(0){}

	bool isvalid()const{
		return xmin <= xmax && ymin <= ymax && zmin <= zmax;
	}
};



template <typename T>
class closest_point_finder : public t_traits<T>{


	std::vector<const T*> points;
	static const unsigned int numdivs = 3;
	std::vector<const T*> pointbins[numdivs+1][numdivs+1][numdivs+1];
	boundingbox<T> bbox;
	real_t binsize;


	struct DistComp{
		const T& point;
		DistComp(const T& point):point(point){}
		bool operator()(const T* lhs, const T* rhs){
			return squared_distance(*lhs, point) < squared_distance(*rhs, point);
		}
	};

	struct DistCompSafe{
		const T& point;
		DistCompSafe(const T& point):point(point){}
		bool operator()(const T*& lhs, const T*& rhs){
			if(lhs && rhs)
				return squared_distance(*lhs, point) < squared_distance(*rhs, point);
			else if(lhs)
				return true;
			else if(rhs)
				return false;
			else 
				return &lhs < &rhs;
		}
	};

	

	inline const T* find_closest_in_range(const T& point, int xmin, int xmax, 
		int ymin, int ymax, int zmin, int zmax)const{
		const T* closepoint = NULL;
		if(xmin <= xmax && ymin <= ymax && zmin <= zmax){
			real_t sqmindist = -1;
			DistComp distcomp(point);
			for (int ii = xmin; ii <= xmax; ++ii){
				for (int jj = ymin; jj <= ymax; ++jj){
					for (int kk = zmin; kk <= zmax; ++kk){
						const std::vector<const T*>& bin = pointbins[ii][jj][kk];
						if( bin.size() > 0 ){
							const T* binclosepoint = *(std::min_element(bin.begin(), bin.end(), distcomp));
							real_t binsqdist = squared_distance(point, *binclosepoint);
							if(sqmindist < 0 || sqmindist < binsqdist){
								closepoint = binclosepoint;
								sqmindist = binsqdist;
							}
						}
					}
				}
			}

		}
		return closepoint;
	}



	const T* find_closest_at_radius(const T& point, unsigned int radius)const{
		int xindx = std::max(std::min((int )((getx(point) - bbox.xmin)/binsize), (int)numdivs), (int )0);
		int yindx = std::max(std::min((int )((gety(point) - bbox.ymin)/binsize), (int)numdivs), (int )0);
		int zindx = std::max(std::min((int )((getz(point) - bbox.zmin)/binsize), (int)numdivs), (int )0);
		
		const int numpoints = 6+12+8;//6sides +12 edges + 8 corners of a cube/cuboid
		const T* closepoints[numpoints];
		int totalboxes = radius > 0 ? numpoints : 1;

		/*int xx[2] = {std::max((int)(xindx - radius), (int)0 ), std::min((int)(xindx + radius), (int)numdivs)};
		int yy[2] = {std::max((int)(yindx - radius), (int)0), std::min((int)(yindx + radius), (int)numdivs)};
		int zz[2] = {std::max((int)(zindx - radius), (int)0), std::min((int)(zindx + radius), (int)numdivs)};

		const int offset[numpoints][6] = {
			//corners	
			{0, 0, 0, 0, 0, 0},	{0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0},	{0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0},	{0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0, 0},	{0, 0, 0, 0, 0, 0},
			//sides
			{1, -1, 1, -1, 0, 0}, {1, -1, 1, -1, 0, 0},
			{1, -1, 0, 0, 1, -1}, {1, -1, 0, 0, 1, -1},
			{0, 0, 1, -1, 1, -1}, {0, 0, 1, -1, 1, -1},
			//edges
			{1, -1, 0, 0, 0, 0}, {1, -1, 0, 0, 0, 0}, {1, -1, 0, 0, 0, 0}, {1, -1, 0, 0, 0, 0},
			{0, 0, 1, -1, 0, 0}, {0, 0, 1, -1, 0, 0}, {0, 0, 1, -1, 0, 0}, {0, 0, 1, -1, 0, 0},
			{0, 0, 0, 0, 1, -1}, {0, 0, 0, 0, 1, -1}, {0, 0, 0, 0, 1, -1}, {0, 0, 0, 0, 1, -1}
		};

		const int index[numpoints][6] = {
			//corners	
			{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 1}, 
			{0, 0, 1, 1, 0, 0}, {0, 0, 1, 1, 1, 1}, 
			{1, 1, 0, 0, 0, 0}, {1, 1, 0, 0, 1, 1}, 
			{1, 1, 1, 1, 0, 0}, {1, 1, 1, 1, 1, 1},
			//sides
			{0, 1, 0, 1, 0, 0},	{0, 1, 0, 1, 1, 1}, 
			{0, 1, 0, 0, 0, 1},	{0, 1, 1, 1, 0, 1}, 
			{0, 0, 0, 1, 0, 1},	{1, 1, 0, 1, 0, 1}, 
			//edges
			{0, 1, 0, 0, 0, 0},	{0, 1, 1, 1, 0, 0}, {0, 1, 0, 0, 1, 1}, {0, 1, 1, 1, 1, 1}, 
			{0, 0, 0, 1, 0, 0}, {1, 1, 0, 1, 0, 0}, {0, 0, 0, 1, 1, 1}, {1, 1, 0, 1, 1, 1}, 
			{0, 0, 0, 0, 0, 1}, {1, 1, 0, 0, 0, 1}, {0, 0, 1, 1, 0, 1}, {1, 1, 1, 1, 0, 1}
		};

		
		for ( int tt = 0; tt < totalboxes; ++tt){

			int xmin = xx[index[tt][0]] + offset[tt][0];
			int xmax = xx[index[tt][1]] + offset[tt][1];
			int ymin = yy[index[tt][2]] + offset[tt][2];
			int ymax = yy[index[tt][3]] + offset[tt][3];
			int zmin = zz[index[tt][4]] + offset[tt][4];
			int zmax = zz[index[tt][5]] + offset[tt][5];
			closepoints[tt] = NULL;
			if(xmin <= xmax && ymin <= ymax && zmin <= zmax){
				real_t sqmindist = -1;
				DistComp distcomp(point);
				for (int ii = xmin; ii <= xmax; ++ii){
					for (int jj = ymin; jj <= ymax; ++jj){
						for (int kk = zmin; kk <= zmax; ++kk){
							const std::vector<const T*>& bin = pointbins[ii][jj][kk];
							if( bin.size() > 0 ){
								const T* binclosepoint = *(std::min_element(bin.begin(), bin.end(), distcomp));
								real_t binsqdist = squared_distance(point, *binclosepoint);
								if(sqmindist < 0 || sqmindist < binsqdist){
									closepoints[tt] = binclosepoint;
									sqmindist = binsqdist;
								}
							}
						}
					}
				}
			}
		}
		*/
		

		int xmin = std::max((int)(xindx - radius), (int)0);
		int xmax = std::min((int)(xindx + radius), (int)numdivs);
		int ymin = std::max((int)(yindx - radius), (int)0);
		int ymax = std::min((int)(yindx + radius), (int)numdivs);
		int zmin = std::max((int)(zindx - radius), (int)0);
		int zmax = std::min((int)(zindx + radius), (int)numdivs);


		if (radius == 0){
			return find_closest_in_range(point, xmin, xmax, ymin, ymax, zmin, zmax);
		}

		closepoints[0] = find_closest_in_range(point, xmin+1, xmax-1, ymin+1, ymax-1, zmin, zmin);
		closepoints[1] = find_closest_in_range(point, xmin+1, xmax-1, ymin+1, ymax-1, zmax, zmax);
		closepoints[2] = find_closest_in_range(point, xmin+1, xmax-1, ymin, ymin, zmin+1, zmax-1);
		closepoints[3] = find_closest_in_range(point, xmin+1, xmax-1, ymax, ymax, zmin+1, zmax-1);
		closepoints[4] = find_closest_in_range(point, xmin, xmin, ymin+1, ymax-1, zmin+1, zmax-1);
		closepoints[5] = find_closest_in_range(point, xmax, xmax, ymin+1, ymax-1, zmin+1, zmax-1);
			
		
		closepoints[6] = find_closest_in_range(point, xmin+1, xmax-1, ymin, ymin, zmin, zmin);
		closepoints[7] = find_closest_in_range(point, xmin+1, xmax-1, ymax, ymax, zmin, zmin);
		closepoints[8] = find_closest_in_range(point, xmin+1, xmax-1, ymin, ymin, zmax, zmax);
		closepoints[9] = find_closest_in_range(point, xmin+1, xmax-1, ymax, ymax, zmax, zmax);
			

		closepoints[10] = find_closest_in_range(point, xmin, xmin, ymin+1, ymax-1, zmin, zmin);
		closepoints[11] = find_closest_in_range(point, xmax, xmax, ymin+1, ymax-1, zmin, zmin);
		closepoints[12] = find_closest_in_range(point, xmin, xmin, ymin+1, ymax-1, zmax, zmax);
		closepoints[13] = find_closest_in_range(point, xmax, xmax, ymin+1, ymax-1, zmax, zmax);

			 

		closepoints[14] = find_closest_in_range(point, xmin, xmin, ymin, ymin, zmin+1, zmax-1);
		closepoints[15] = find_closest_in_range(point, xmax, xmax, ymin, ymin, zmin+1, zmax-1);
		closepoints[16] = find_closest_in_range(point, xmin, xmin, ymax, ymax, zmin+1, zmax-1);
		closepoints[17] = find_closest_in_range(point, xmax, xmax, ymax, ymax, zmin+1, zmax-1);
		
		
		closepoints[18] = find_closest_in_range(point, xmin, xmin, ymin, ymin, zmin, zmin);
		closepoints[19] = find_closest_in_range(point, xmin, xmin, ymin, ymin, zmax, zmax);
		closepoints[20] = find_closest_in_range(point, xmin, xmin, ymax, ymax, zmin, zmin);
		closepoints[21] = find_closest_in_range(point, xmin, xmin, ymax, ymax, zmax, zmax);
		closepoints[22] = find_closest_in_range(point, xmax, xmax, ymin, ymin, zmin, zmin);
		closepoints[23] = find_closest_in_range(point, xmax, xmax, ymin, ymin, zmax, zmax);
		closepoints[24] = find_closest_in_range(point, xmax, xmax, ymax, ymax, zmin, zmin);
		closepoints[25] = find_closest_in_range(point, xmax, xmax, ymax, ymax, zmax, zmax);

		DistCompSafe distcomp(point);
		const T* closepoint = *(std::min_element(closepoints, closepoints + totalboxes, distcomp));
		return closepoint;

	}


public:

	void recalculate(){
		bbox = boundingbox<T>(points.begin(), points.end());
		
		binsize = std::max(std::max((bbox.xmax - bbox.xmin)/numdivs, (bbox.ymax - bbox.ymin)/numdivs), (bbox.zmax - bbox.zmin)/numdivs);
		for(std::vector<const T*>::const_iterator iter = points.begin(); iter != points.end() ; ++iter){
			unsigned int xindx = (unsigned int )((getx(**iter) - bbox.xmin)/binsize);
			unsigned int yindx = (unsigned int )((gety(**iter) - bbox.ymin)/binsize);
			unsigned int zindx = (unsigned int )((getz(**iter) - bbox.zmin)/binsize);
			pointbins[xindx][yindx][zindx].push_back(*iter);
		}
	}



	
	
	const T* find_closest(const T& point)const{

		if ( points.size() == 0 )
			return NULL;
		bool debugit = false;
		//DWORD t1, t2, t3;
		//t1 = GetTickCount();

		unsigned int radius = 0;
		const T* closepoint = NULL;
		while(!closepoint) {
			//loop in increasing radii until we find a close point.
			closepoint = find_closest_at_radius(point, radius++);
		}

		real_t sqmindist1 = squared_distance(*closepoint, point);

		int xindx = std::max(std::min((int )((getx(point) - bbox.xmin)/binsize), (int)numdivs), (int )0);
		int yindx = std::max(std::min((int )((gety(point) - bbox.ymin)/binsize), (int)numdivs), (int )0);
		int zindx = std::max(std::min((int )((getz(point) - bbox.zmin)/binsize), (int)numdivs), (int )0);

		
		T centerpoint = t_traits<T>::makeT<T>( bbox.xmin + xindx*binsize + binsize/2, 
						bbox.ymin + yindx*binsize + binsize/2, bbox.zmin + zindx*binsize + binsize/2 );
		real_t sqmindist2 = squared_distance(*closepoint, centerpoint);
		real_t maxradiusf = sqrt(std::max(sqmindist1, sqmindist2))/binsize;

		unsigned int maxradius = (unsigned int)(maxradiusf) + (maxradiusf > floor(maxradiusf) ? 1 : 0);
		while (  radius < maxradius ){

			const T* another_closepoint = find_closest_at_radius(point, radius++);
			if (squared_distance(*another_closepoint, point) < sqmindist1){
				closepoint = another_closepoint;
				sqmindist1 = squared_distance(*closepoint, point);
				sqmindist2 = squared_distance(*closepoint, centerpoint);

				real_t maxradiusf = sqrt(std::max(sqmindist1, sqmindist2))/binsize;
				maxradius = (unsigned int)(maxradiusf) + (maxradiusf > floor(maxradiusf) ? 1 : 0);
			}
			//debugit = true;
		}
		
		//t2 =  GetTickCount();
		//if(debugit)
			//cout << "POINT " << point << " TIME = " << t2-t1 << endl;
		return closepoint;
	}

	const T* find_closest_linear(const T& point)const{
		const T* closestpoint = NULL;
		if ( points.size() > 0 ){
			real_t mindist =  distance(*points.front(), point);
			closestpoint = points.front();
			for(std::vector<const T*>::const_iterator piter = points.begin() + 1; piter != points.end(); ++piter){
				real_t curdist = distance(**piter, point);
				closestpoint = curdist < mindist ? *piter : closestpoint;
				mindist = std::min(mindist, curdist);
			}
		}
		return closestpoint;
	}

	//input is begin and end iterators to a container of pointer of Ts
	template <typename IterT> 		
	closest_point_finder(IterT beginP, IterT endP,
		typename std::enable_if< std::is_same<typename std::iterator_traits<IterT>::value_type, T* >::value >::type* dummy = 0)
	{
		for(IterT it = beginP; it != endP; ++it){
			points.push_back((*it));
		}
		recalculate();
	}

	//input is begin and end iterators to a container of const pointer of Ts
	template <typename IterT> 		
	closest_point_finder(IterT beginP, IterT endP,
		typename std::enable_if< std::is_same<typename std::iterator_traits<IterT>::value_type, const T* >::value >::type* dummy = 0)
	{
		for(IterT it = beginP; it != endP; ++it){
			points.push_back((*it));
		}
		recalculate();
	}

	//input is begin and end iterators to a container of Ts
	template <typename IterT> 		
	closest_point_finder(IterT beginP, IterT endP,
		typename std::enable_if< std::is_same<typename std::iterator_traits<IterT>::value_type, T >::value >::type* dummy = 0)		
	{	
		for(IterT it = beginP; it != endP; ++it){
			points.push_back(&(*it));
		}
		recalculate();
	}

	//input is a 
	template <typename IterT> 		
	closest_point_finder(IterT beginP, IterT endP,
		typename std::enable_if< std::is_same<typename std::iterator_traits<IterT>::value_type::first_type, const T >::value >::type* dummy = 0)		
	{	
		for(IterT it = beginP; it != endP; ++it){
			points.push_back(&(it->first));
		}
		recalculate();
	}

	template <typename IterT> 		
	closest_point_finder(IterT beginP, IterT endP,
		typename std::enable_if< std::is_same<typename std::iterator_traits<IterT>::value_type::second_type, T >::value >::type* dummy = 0)		
	{	
		for(IterT it = beginP; it != endP; ++it){
			points.push_back(&(it->second));
		}
		recalculate();
	}

	closest_point_finder(const T* start, const T* end)
	{
		for(const T* iter = start; iter != end; ++iter){
			points.push_back(iter);
		}
		recalculate();
	}

	closest_point_finder(){}
};//SpaceHash


}//namespace

#endif


/* 
	typedef float (T::*FloatType)()const;
	typedef double (T::*DoubleType)()const;


	template <typename std::enable_if< has_x<T, DoubleType>::value>::type>
	struct realdef{
		typedef double real_t;
	};

	template <typename std::enable_if< has_x<T, FloatType>::value>::type>
	struct realdef{
		typedef float real_t;
	};

	typedef realdef::real_t real_t;*/
