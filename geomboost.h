
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

	typedef T* TRef;

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
};



template <typename T>
class closest_point_finder : public t_traits<T>{

	std::vector<const T*> points;
	static const unsigned int numdivs = 3;
	std::vector<const T*> binpoints[numdivs+1][numdivs+1][numdivs+1];
	boundingbox<T> bbox;
	real_t binx;
	real_t biny;
	real_t binz;


	struct DistComp{
		const T& point;
		DistComp(const T& point):point(point){}
		bool operator()(const T* lhs, const T* rhs){
			return squared_distance(*lhs, point) < squared_distance(*rhs, point);
		}
	};


public:

	void recalculate(){
		bbox = boundingbox<T>(points.begin(), points.end());
		binx = (bbox.xmax - bbox.xmin)/numdivs;
		biny = (bbox.ymax - bbox.ymin)/numdivs;
		binz = (bbox.zmax - bbox.zmin)/numdivs;
		for(std::vector<const T*>::const_iterator iter = points.begin(); iter != points.end() ; ++iter){
			unsigned int xindx = (unsigned int )((getx(**iter) - bbox.xmin)/binx);
			unsigned int yindx = (unsigned int )((gety(**iter) - bbox.ymin)/biny);
			unsigned int zindx = (unsigned int )((getz(**iter) - bbox.zmin)/binz);
			binpoints[xindx][yindx][zindx].push_back(*iter);
		}
	}

	
	const T* find_closest(const T& point)const{
		unsigned int xindx = std::max(std::min((unsigned int )((getx(point) - bbox.xmin)/binx), numdivs), (unsigned int )0);
		unsigned int yindx = std::max(std::min((unsigned int )((gety(point) - bbox.ymin)/biny), numdivs), (unsigned int )0);
		unsigned int zindx = std::max(std::min((unsigned int )((getz(point) - bbox.zmin)/binz), numdivs), (unsigned int )0);
		
		DistComp distcomp(point);
		const std::vector<const T*>& bin = binpoints[xindx][yindx][zindx];
		return bin.size() == 0 ? NULL : *(std::min_element(bin.begin(), bin.end(), distcomp));
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
