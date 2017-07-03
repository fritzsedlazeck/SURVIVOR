################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/vcfs/Annotate_vcf.cpp \
../src/vcfs/Combine_3_VCF.cpp \
../src/vcfs/Compoverlap_VCF.cpp \
../src/vcfs/Detect_nested.cpp \
../src/vcfs/Filter_vcf.cpp \
../src/vcfs/Merge_VCF.cpp 

OBJS += \
./src/vcfs/Annotate_vcf.o \
./src/vcfs/Combine_3_VCF.o \
./src/vcfs/Compoverlap_VCF.o \
./src/vcfs/Detect_nested.o \
./src/vcfs/Filter_vcf.o \
./src/vcfs/Merge_VCF.o 

CPP_DEPS += \
./src/vcfs/Annotate_vcf.d \
./src/vcfs/Combine_3_VCF.d \
./src/vcfs/Compoverlap_VCF.d \
./src/vcfs/Detect_nested.d \
./src/vcfs/Filter_vcf.d \
./src/vcfs/Merge_VCF.d 


# Each subdirectory must supply rules for building sources it contributes
src/vcfs/%.o: ../src/vcfs/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


