################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/snp_overlap/Overlap_snps.cpp 

OBJS += \
./src/snp_overlap/Overlap_snps.o 

CPP_DEPS += \
./src/snp_overlap/Overlap_snps.d 


# Each subdirectory must supply rules for building sources it contributes
src/snp_overlap/%.o: ../src/snp_overlap/%.cpp src/snp_overlap/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


