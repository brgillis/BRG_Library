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

class my_fftw_plan
{
private:
	fftw_plan _p_;
public:
	// Default constructor
	my_fftw_plan() : _p_(nullptr) {}

	// Construct from fftw_plan
	my_fftw_plan(const fftw_plan & p) : _p_(p) {}
	my_fftw_plan(fftw_plan && p) : _p_(std::move(p)) {}

	// Destructor/destroyer
	~my_fftw_plan() { fftw_destroy_plan(_p_); }
	void destroy()
	{
		fftw_destroy_plan(_p_);
		_p_ = nullptr;
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
		fftw_import_wisdom_from_filename(_filename_.c_str());
	}
	fftw_wisdom_accumulator(const std::string & filename)
	: fftw_wisdom_accumulator(std::string(filename))
	{
	}
	fftw_wisdom_accumulator()
	: fftw_wisdom_accumulator(std::string(".fftw_wisdom"))
	{
	}

	void load(const std::string & filename)
	{
		fftw_import_wisdom_from_filename(filename.c_str());
	}

	void save(const std::string & filename)
	{
		fftw_export_wisdom_to_filename(filename.c_str());
	}

	~fftw_wisdom_accumulator()
	{
		fftw_export_wisdom_to_filename(_filename_.c_str());
	}
};

// Typedefs

#if(1)

typedef double flt_type;
typedef Eigen::Array<flt_type,1,Eigen::Dynamic> flt_array_type;
typedef Eigen::Map<flt_array_type> flt_array_map_type;

typedef std::complex<flt_type> complex_type;
typedef Eigen::Array<complex_type,1,Eigen::Dynamic> complex_array_type;
typedef Eigen::Map<complex_array_type> complex_array_map_type;

typedef std::unique_ptr<flt_type,decltype(fftw_free)> fftw_flt_ptr;
typedef std::unique_ptr<fftw_complex,decltype(fftw_free)> fftw_complex_ptr;

typedef std::unique_ptr<fftw_plan,decltype(fftw_destroy_plan)> fftw_plan_ptr;

#endif

// Fourier transform of discrete values
#if(1)
template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
complex_array_type Fourier_transform(array_type && vals)
{
	int N = ssize(vals);
	fftw_complex_ptr out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	#pragma openmp critical(make_fftw_plan)
	my_fftw_plan plan = fftw_plan_dft_r2c_1d(N, vals.data(), out, FFTW_ESTIMATE);

	plan.execute();

	complex_array_type result(N);

	for(int i=0; i<N; ++i)
	{
		result[i] = complex_type(out[i][0],out[i][1]);
	}

	return result;
}
template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
complex_array_type Fourier_transform(const array_type & vals)
{
	return Fourier_transform(array_type(vals));
}

template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
complex_array_type planned_Fourier_transform(array_type && vals,
		boost::optional<my_fftw_plan> & plan)
{
	int N = ssize(vals);
	fftw_complex_ptr out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	// If we don't have a plan set up, make one now
	if(!plan)
	{
		// Use a temporary object so the values won't be overwritten when the plan is
		// measured
		#pragma openmp critical(make_fftw_plan)
		{
			fftw_wisdom_accumulator wisdom_accumulator;
			plan = fftw_plan_dft_r2c_1d(N, array_type(N).data(), out.get(), FFTW_MEASURE);
		}
	}

	// Execute the plan on the proper data
	fftw_execute_dft_r2c(plan->get_plan(),vals.data(),out.get());

	complex_array_type result(N);

	for(int i=0; i<N; ++i)
	{
		result[i] = complex_type(out[i][0],out[i][1]);
	}

	return result;
}
template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
complex_array_type planned_Fourier_transform(const array_type & vals,
		boost::optional<my_fftw_plan> & plan)
{
	return planned_Fourier_transform(array_type(vals), plan);
}

#endif // Fourier transform of discrete values

// Inverse Fourier transform of discrete values
#if(1)
template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
flt_array_type inverse_Fourier_transform(array_type && vals)
{
	int N = ssize(vals);

	fftw_complex_ptr in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	for(int i=0; i<N; ++i)
	{
		in[i][0] = vals[i].real();
		in[i][1] = vals[i].imag();
	}

	fftw_flt_ptr out = (flt_type*) fftw_malloc(sizeof(flt_type) * N);

	#pragma openmp critical(make_fftw_plan)
	my_fftw_plan plan = fftw_plan_dft_c2r_1d(N, in.get(), out.get(), FFTW_ESTIMATE);

	plan.execute();

	flt_array_type result = flt_array_map_type(N,out.get())/N;

	return result;
}
template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
flt_array_type inverse_Fourier_transform(const array_type & vals)
{
	return Fourier_inverse_transform(array_type(vals));
}

template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
flt_array_type planned_inverse_Fourier_transform(array_type && vals,
		boost::optional<my_fftw_plan> & plan)
{
	int N = ssize(vals);

	fftw_complex_ptr in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	fftw_flt_ptr out = (flt_type*) fftw_malloc(sizeof(flt_type) * N);

	// If we don't have a plan set up, make one now
	if(!plan)
	{
		#pragma openmp critical(make_fftw_plan)
		{
			fftw_wisdom_accumulator wisdom_accumulator;
			plan = fftw_plan_dft_c2r_1d(N, in.get(), out.get(), FFTW_MEASURE);
		}
	}

	// Now fill in the array (this must be after the plan, as the plan may overwrite it)
	for(int i=0; i<N; ++i)
	{
		in[i][0] = vals[i].real();
		in[i][1] = vals[i].imag();
	}

	// Execute the plan on the proper data
	fftw_execute_dft_c2r(plan->get_plan(),in.get(),out.get());

	flt_array_type result = flt_array_map_type(N,out.get())/N;

	return result;
}
template< typename array_type,
typename std::enable_if<brgastro::is_stl_or_eigen_container<array_type>::value,char>::type = 0 >
flt_array_type planned_inverse_Fourier_transform(const array_type & vals,
		boost::optional<my_fftw_plan> & plan)
{
	return planned_inverse_Fourier_transform(array_type(vals), plan);
}

#endif // Inverse Fourier transform of discrete values

// Fourier transform of function
#if(1)

template< typename func_type,
typename std::enable_if<!brgastro::is_stl_or_eigen_container<func_type>::value,char>::type = 0 >
complex_array_type Fourier_transform( const func_type & func,
		const flt_type & min = 0.,
		const flt_type & max = 1.,
		const int & samples = 1000)
{
	flt_array_type vals = flt_array_type::LinSpaced(samples,min,max).UnaryExpr(func);
	return Fourier_transform(std::move(vals));
}

template< typename func_type,
typename std::enable_if<!brgastro::is_stl_or_eigen_container<func_type>::value,char>::type = 0 >
complex_array_type planned_Fourier_transform( const func_type & func,
		boost::optional<my_fftw_plan> plan,
		const flt_type & min = 0.,
		const flt_type & max = 1.,
		const int & samples = 1000)
{
	flt_array_type vals = flt_array_type::LinSpaced(samples,min,max).UnaryExpr(func);
	return planned_Fourier_transform(std::move(vals),plan);
}


#endif

} // namespace Fourier

} // namespace brgastro

#endif /* BRG_MATH_FOURIER_HPP_ */
