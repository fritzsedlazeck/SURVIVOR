################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/scrup_vcf/scrup_svs.cpp 

OBJS += \
./src/scrup_vcf/scrup_svs.o 

CPP_DEPS += \
./src/scrup_vcf/scrup_svs.d 


# Each subdirectory must supply rules for building sources it contributes
src/scrup_vcf/%.o: ../src/scrup_vcf/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


