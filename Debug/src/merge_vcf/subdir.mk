################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/merge_vcf/IntervallTree.cpp \
../src/merge_vcf/combine_svs.cpp 

OBJS += \
./src/merge_vcf/IntervallTree.o \
./src/merge_vcf/combine_svs.o 

CPP_DEPS += \
./src/merge_vcf/IntervallTree.d \
./src/merge_vcf/combine_svs.d 


# Each subdirectory must supply rules for building sources it contributes
src/merge_vcf/%.o: ../src/merge_vcf/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


