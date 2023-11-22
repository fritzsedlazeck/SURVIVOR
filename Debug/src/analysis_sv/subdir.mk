################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/analysis_sv/Density_VCF.cpp \
../src/analysis_sv/GIAB_summary.cpp \
../src/analysis_sv/MT_identifier.cpp \
../src/analysis_sv/MUMmer_overlap.cpp \
../src/analysis_sv/Select_samples.cpp \
../src/analysis_sv/Simplify_SVs.cpp \
../src/analysis_sv/Summ_mat.cpp 

OBJS += \
./src/analysis_sv/Density_VCF.o \
./src/analysis_sv/GIAB_summary.o \
./src/analysis_sv/MT_identifier.o \
./src/analysis_sv/MUMmer_overlap.o \
./src/analysis_sv/Select_samples.o \
./src/analysis_sv/Simplify_SVs.o \
./src/analysis_sv/Summ_mat.o 

CPP_DEPS += \
./src/analysis_sv/Density_VCF.d \
./src/analysis_sv/GIAB_summary.d \
./src/analysis_sv/MT_identifier.d \
./src/analysis_sv/MUMmer_overlap.d \
./src/analysis_sv/Select_samples.d \
./src/analysis_sv/Simplify_SVs.d \
./src/analysis_sv/Summ_mat.d 


# Each subdirectory must supply rules for building sources it contributes
src/analysis_sv/%.o: ../src/analysis_sv/%.cpp src/analysis_sv/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


