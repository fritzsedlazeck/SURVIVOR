################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/simulator/Eval_vcf.cpp \
../src/simulator/Pac_Simulator.cpp \
../src/simulator/SV_Simulator.cpp 

OBJS += \
./src/simulator/Eval_vcf.o \
./src/simulator/Pac_Simulator.o \
./src/simulator/SV_Simulator.o 

CPP_DEPS += \
./src/simulator/Eval_vcf.d \
./src/simulator/Pac_Simulator.d \
./src/simulator/SV_Simulator.d 


# Each subdirectory must supply rules for building sources it contributes
src/simulator/%.o: ../src/simulator/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


