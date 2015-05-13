/*
 * Fourier.hpp
 *
 *  Created on: 11 May 2015
 *      Author: brg
 */

#ifndef BRG_MATH_FOURIER_HPP_
#define BRG_MATH_FOURIER_HPP_

#include <complex>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <boost/optional.hpp>
#include <fftw3.h>
#include <Eigen/Core>

#include "brg/container/is_eigen_container.hpp"
#include "brg/utility.hpp"

namespace brgastro {

namespace Fourier {

constexpr const char * default_fftw_wisdom_filename = ".fftw_wisdom";

class my_fftw_plan
{
private:
	fftw_plan _p_;
public:
	// Default constructor
	my_fftw_plan() : _p_(nullptr) {}

	// Construct from fftw_plan
	my_fftw_plan(const fftw_plan & p) : _p_(p) {}

	// Destructor/destroyer
	~my_fftw_plan() { fftw_destroy_plan(_p_); }
	void destroy()
	{
		fftw_destroy_plan(_p_);
		_p_ = nullptr;
	}

	// Operator=
	my_fftw_plan & operator=(const my_fftw_plan &) = delete;
	my_fftw_plan & operator=(my_fftw_plan && other)
	{
		std::swap(_p_,other._p_);
		return *this;
	}

	// Execute
	void execute() {return fftw_execute(_p_);}

	// Convert to fftw_plan
	fftw_plan & get_plan() noexcept {return _p_;}
	operator fftw_plan &() noexcept {return _p_;}
};

class fftw_wisdom_accumulator
{
private:
	std::string _filename_;
public:
	fftw_wisdom_accumulator(std::string && filename)
	: _filename_(std::move(filename))
	{
		// Load wisdom unless we have a null filename
		if(_filename_.size()>0)
		{
			fftw_import_wisdom_from_filename(_filename_.c_str()); // Ignore failure
		}
	}
	fftw_wisdom_accumulator(const std::string & filename)
	: fftw_wisdom_accumulator(std::string(filename))
	{
	}
	fftw_wisdom_accumulator()
	: fftw_wisdom_accumulator(std::string(default_fftw_wisdom_filename))
	{
	}

	void load(const std::string & filename)
	{
		short res;

		#pragma omp critical(load_or_save_fftw_wisdom)
		res = fftw_import_wisdom_from_filename(filename.c_str());

		// Throw an exception if we couldn't load
		if(!res)
		{
			throw std::runtime_error("Could not load wisdom from " + filename + ".");
		}
	}

	void save(const std::string & filename)
	{
		short res;

		#pragma omp critical(load_or_save_fftw_wisdom)
		res = fftw_export_wisdom_to_filename(filename.c_str());

		// Throw an exception if we couldn't load
		if(!res)
		{
			throw std::runtime_error("Could not save wisdom to " + filename + ".");
		}
	}

	~fftw_wisdom_accumulator()
	{
		// When this goes out of scope, save it unless it has a null filename
		if(_filename_.size()>0)
		{
			#pragma omp critical(load_or_save_fftw_wisdom)
			fftw_export_wisdom_to_filename(_filename_.c_str()); // Ignore failure
		}
	}
};

template< typename T >
struct fftw_array_deleter
{
	void operator()(T* p)
	{
		fftw_free(p);
	}
};

// Typedefs

#if(1)

typedef double flt_type;
typedef Eigen::Array<flt_type,Eigen::Dynamic,1> flt_array_type;
typedef Eigen::Map<flt_array_type> flt_array_map_type;

typedef std::complex<flt_type> complex_type;
typedef Eigen::Array<complex_type,Eigen::Dynamic,1> complex_array_type;
typedef Eigen::Map<complex_array_type> complex_array_map_type;

typedef std::unique_ptr<flt_type,fftw_array_deleter<flt_type>> fftw_flt_ptr;
typedef std::unique_ptr<fftw_complex,fftw_array_deleter<fftw_complex>> fftw_complex_ptr;

#endif

// Fourier transform of discrete values
#if(1)
template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
complex_array_type Fourier_transform(const array_type & vals,
		const boost::optional<fftw_wisdom_accumulator> & wisdom = boost::none )
{
	int N = ssize(vals);
	fftw_complex_ptr out((fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N));

	my_fftw_plan plan;

	// Check if we have any wisdom. If not, estimate. If so, measure.
	if(!wisdom)
	{
		// For the estimating case, it's safe to use the vals array as input - it won't be overwritten.
		// We do have to use a const_cast so it can go into the function though
		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_dft_r2c_1d(N, const_cast<flt_type *>(vals.data()), out.get(), FFTW_ESTIMATE);

		plan.execute();
	}
	else
	{
		fftw_flt_ptr in((flt_type*) fftw_malloc(sizeof(flt_type) * N));

		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_dft_r2c_1d(N, in.get(), out.get(), FFTW_MEASURE);

		for(int i=0; i<N; ++i) in.get()[i] = vals[i]; // Need to initialize after planning in this case

		plan.execute(); // Execute before in is allocated
	}

	complex_array_type result(N);

	for(int i=0; i<N; ++i)
	{
		result(i) = complex_type(out.get()[i][0],out.get()[i][1]);
	}

	return result;
}

#endif // Fourier transform of discrete values

// Inverse Fourier transform of discrete values
#if(1)
template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
flt_array_type inverse_Fourier_transform(const array_type & vals,
		const boost::optional<fftw_wisdom_accumulator> & wisdom = boost::none )
{
	int N = ssize(vals);

	fftw_complex_ptr in((fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N));
	fftw_flt_ptr out((flt_type*) fftw_malloc(sizeof(flt_type) * N));

	my_fftw_plan plan;

	// Check if we have any wisdom. If not, estimate. If so, measure.
	if(!wisdom)
	{
		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_dft_c2r_1d(N, in.get(), out.get(), FFTW_ESTIMATE);
	}
	else
	{
		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_dft_c2r_1d(N, in.get(), out.get(), FFTW_MEASURE);
	}

	// Initialize after planning
	for(int i=0; i<N; ++i)
	{
		in.get()[i][0] = vals[i].real();
		in.get()[i][1] = vals[i].imag();
	}

	plan.execute();

	flt_array_type result = flt_array_map_type(out.get(),N)/N;

	return result;
}

#endif // Inverse Fourier transform of discrete values

// Fourier sin transform of discrete values
#if(1)
template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
flt_array_type Fourier_sin_transform(const array_type & vals,
		const boost::optional<fftw_wisdom_accumulator> & wisdom = boost::none )
{
	int N = ssize(vals);
	fftw_flt_ptr out((flt_type*) fftw_malloc(sizeof(flt_type) * N));

	my_fftw_plan plan;

	// Check if we have any wisdom. If not, estimate. If so, measure.
	if(!wisdom)
	{
		// For the estimating case, it's safe to use the vals array as input - it won't be overwritten.
		// We do have to use a const_cast so it can go into the function though
		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_r2r_1d(N, const_cast<flt_type *>(vals.data()), out.get(),
				FFTW_RODFT10, FFTW_ESTIMATE);

		plan.execute();
	}
	else
	{
		fftw_flt_ptr in((flt_type*) fftw_malloc(sizeof(flt_type) * N));

		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_r2r_1d(N, in.get(), out.get(), FFTW_RODFT10, FFTW_MEASURE);

		for(int i=0; i<N; ++i) in.get()[i] = vals[i]; // Need to initialize after planning in this case

		plan.execute(); // Execute before in is deallocated
	}

	flt_array_type result = flt_array_map_type(out.get(),N);

	return result;
}

#endif // Fourier sin transform of discrete values

// Inverse Fourier sin transform of discrete values
#if(1)

template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
flt_array_type inverse_Fourier_sin_transform(const array_type & vals,
		const boost::optional<fftw_wisdom_accumulator> & wisdom = boost::none )
{
	int N = ssize(vals);

	fftw_flt_ptr out((flt_type*) fftw_malloc(sizeof(flt_type) * N));

	my_fftw_plan plan;

	// Check if we have any wisdom. If not, estimate. If so, measure.
	if(!wisdom)
	{
		// For the estimating case, it's safe to use the vals array as input - it won't be overwritten.
		// We do have to use a const_cast so it can go into the function though
		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_r2r_1d(N, const_cast<flt_type *>(vals.data()), out.get(),
				FFTW_RODFT01, FFTW_ESTIMATE);

		plan.execute();
	}
	else
	{
		fftw_flt_ptr in((flt_type*) fftw_malloc(sizeof(flt_type) * N));

		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_r2r_1d(N, in.get(), out.get(), FFTW_RODFT01, FFTW_MEASURE);

		for(int i=0; i<N; ++i) in.get()[i] = vals[i]; // Need to initialize after planning in this case

		plan.execute(); // Execute before in is deallocated
	}

	flt_array_type result = flt_array_map_type(out.get(),N)/(2*N);

	return result;
}

#endif // Inverse Fourier sin transform of discrete values

// spherical Fourier transform of discrete values
#if(1)
template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
flt_array_type spherical_Fourier_transform(array_type vals,
		const boost::optional<fftw_wisdom_accumulator> & wisdom = boost::none )
{
	int N = ssize(vals);
	fftw_flt_ptr out((flt_type*) fftw_malloc(sizeof(flt_type) * N));

	vals *= flt_array_type::LinSpaced(N,0.5,N-0.5);

	my_fftw_plan plan;

	// Check if we have any wisdom. If not, estimate. If so, measure.
	if(!wisdom)
	{
		// For the estimating case, it's safe to use the vals array as input - it won't be overwritten.
		// We do have to use a const_cast so it can go into the function though
		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_r2r_1d(N, const_cast<flt_type *>(vals.data()), out.get(),
				FFTW_RODFT10, FFTW_ESTIMATE);

		plan.execute();
	}
	else
	{
		fftw_flt_ptr in((flt_type*) fftw_malloc(sizeof(flt_type) * N));

		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_r2r_1d(N, in.get(), out.get(), FFTW_RODFT10, FFTW_MEASURE);

		for(int i=0; i<N; ++i) in.get()[i] = vals[i]; // Need to initialize after planning in this case

		plan.execute(); // Execute before in is deallocated
	}

	flt_array_type result = flt_array_map_type(out.get(),N)/flt_array_type::LinSpaced(N,1,N);

	return result;
}

#endif // spherical Fourier transform of discrete values

// Inverse spherical Fourier transform of discrete values
#if(1)
template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
flt_array_type inverse_spherical_Fourier_transform(array_type vals,
		const boost::optional<fftw_wisdom_accumulator> & wisdom = boost::none )
{
	int N = ssize(vals);

	fftw_flt_ptr in((flt_type*) fftw_malloc(sizeof(flt_type) * N));
	fftw_flt_ptr out((flt_type*) fftw_malloc(sizeof(flt_type) * N));

	vals *= flt_array_type::LinSpaced(N,1,N);

	my_fftw_plan plan;

	// Check if we have any wisdom. If not, estimate. If so, measure.
	if(!wisdom)
	{
		// For the estimating case, it's safe to use the vals array as input - it won't be overwritten.
		// We do have to use a const_cast so it can go into the function though
		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_r2r_1d(N, const_cast<flt_type *>(vals.data()), out.get(),
				FFTW_RODFT01, FFTW_ESTIMATE);

		plan.execute();
	}
	else
	{
		fftw_flt_ptr in((flt_type*) fftw_malloc(sizeof(flt_type) * N));

		#pragma omp critical(make_fftw_plan)
		plan = fftw_plan_r2r_1d(N, in.get(), out.get(), FFTW_RODFT01, FFTW_MEASURE);

		for(int i=0; i<N; ++i) in.get()[i] = vals[i]; // Need to initialize after planning in this case

		plan.execute(); // Execute before in is deallocated
	}

	flt_array_type result = flt_array_map_type(out.get(),N)/(2*N*flt_array_type::LinSpaced(N,0.5,N-0.5));

	return result;
}

#endif // Inverse spherical Fourier transform of discrete values

// Fourier transform of function
#if(1)

template< typename func_type,
typename std::enable_if<!brgastro::is_stl_or_eigen_container<func_type>::value,char>::type = 0 >
complex_array_type Fourier_transform( const func_type & func,
		const flt_type & min = 0.,
		const flt_type & max = 1.,
		const int & samples = 1024,
		const boost::optional<fftw_wisdom_accumulator> & wisdom = boost::none)
{
	flt_array_type vals = flt_array_type::LinSpaced(samples,min,max).unaryExpr(func);
	return Fourier_transform(std::move(vals),wisdom);
}

template< typename func_type,
typename std::enable_if<!brgastro::is_stl_or_eigen_container<func_type>::value,char>::type = 0 >
flt_array_type Fourier_sin_transform( const func_type & func,
		const flt_type & max = 1.,
		const int & samples = 1024,
		const boost::optional<fftw_wisdom_accumulator> & wisdom = boost::none)
{
	flt_array_type xvals = flt_array_type::LinSpaced(samples,0.,max) + max/(2*(samples-1));
	flt_array_type vals = xvals.unaryExpr(func);

	return Fourier_sin_transform(std::move(vals),wisdom);
}

template< typename func_type,
typename std::enable_if<!brgastro::is_stl_or_eigen_container<func_type>::value,char>::type = 0 >
flt_array_type spherical_Fourier_transform( const func_type & func,
		const flt_type & max = 1.,
		const int & samples = 1024,
		const boost::optional<fftw_wisdom_accumulator> & wisdom = boost::none)
{
	flt_array_type xvals = flt_array_type::LinSpaced(samples,0.,max) + max/(2*(samples-1));
	flt_array_type vals = xvals * xvals.unaryExpr(func);

	flt_array_type result = (samples-1)/max * Fourier_sin_transform(std::move(vals),wisdom) /
			flt_array_type::LinSpaced(samples,1,samples);

	return result;
}

#endif

} // namespace Fourier

using Fourier::Fourier_transform;
using Fourier::inverse_Fourier_transform;
using Fourier::Fourier_sin_transform;
using Fourier::inverse_Fourier_sin_transform;
using Fourier::spherical_Fourier_transform;
using Fourier::inverse_spherical_Fourier_transform;
using Fourier::fftw_wisdom_accumulator;

} // namespace brgastro

#endif /* BRG_MATH_FOURIER_HPP_ */
