################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../brg_stripping/gabdt.cpp \
../brg_stripping/solve_rt_functors.cpp \
../brg_stripping/stripping_orbit.cpp \
../brg_stripping/stripping_orbit_segment.cpp 

OBJS += \
./brg_stripping/gabdt.o \
./brg_stripping/solve_rt_functors.o \
./brg_stripping/stripping_orbit.o \
./brg_stripping/stripping_orbit_segment.o 

CPP_DEPS += \
./brg_stripping/gabdt.d \
./brg_stripping/solve_rt_functors.d \
./brg_stripping/stripping_orbit.d \
./brg_stripping/stripping_orbit_segment.d 


# Each subdirectory must supply rules for building sources it contributes
brg_stripping/%.o: ../brg_stripping/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -D__GXX_EXPERIMENTAL_CXX0X__ -DNDEBUG -I"/disk2/brg/git/brg_library/brg" -I/disk2/brg/include -O3 -Wall -c -fmessage-length=0 -fopenmp -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


