################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/phasing/Phasing_vcf.cpp 

OBJS += \
./src/phasing/Phasing_vcf.o 

CPP_DEPS += \
./src/phasing/Phasing_vcf.d 


# Each subdirectory must supply rules for building sources it contributes

src/phasing/%.o: ../src/phasing/%.cpp src/phasing/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


