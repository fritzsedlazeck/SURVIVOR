################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CorrectAllele.cpp \
../src/DetectDif.cpp \
../src/Extract_Seq.cpp \
../src/SURVIVOR.cpp \
../src/Summarize_SV.cpp 

OBJS += \
./src/CorrectAllele.o \
./src/DetectDif.o \
./src/Extract_Seq.o \
./src/SURVIVOR.o \
./src/Summarize_SV.o 

CPP_DEPS += \
./src/CorrectAllele.d \
./src/DetectDif.d \
./src/Extract_Seq.d \
./src/SURVIVOR.d \
./src/Summarize_SV.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


