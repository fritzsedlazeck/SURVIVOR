################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/simulator/Error_scanner.cpp \
../src/simulator/Eval_vcf.cpp \
../src/simulator/Pac_Simulator.cpp \
../src/simulator/SV_Simulator.cpp \
../src/simulator/Sim_reads.cpp \
../src/simulator/test_cov.cpp 

OBJS += \
./src/simulator/Error_scanner.o \
./src/simulator/Eval_vcf.o \
./src/simulator/Pac_Simulator.o \
./src/simulator/SV_Simulator.o \
./src/simulator/Sim_reads.o \
./src/simulator/test_cov.o 

CPP_DEPS += \
./src/simulator/Error_scanner.d \
./src/simulator/Eval_vcf.d \
./src/simulator/Pac_Simulator.d \
./src/simulator/SV_Simulator.d \
./src/simulator/Sim_reads.d \
./src/simulator/test_cov.d 


# Each subdirectory must supply rules for building sources it contributes
src/simulator/%.o: ../src/simulator/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


