################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/convert/ConvertMQ0Bed.cpp \
../src/convert/Convert_Pindel.cpp \
../src/convert/Convert_VCF_to_BED.cpp \
../src/convert/Process_Lumpy.cpp 

OBJS += \
./src/convert/ConvertMQ0Bed.o \
./src/convert/Convert_Pindel.o \
./src/convert/Convert_VCF_to_BED.o \
./src/convert/Process_Lumpy.o 

CPP_DEPS += \
./src/convert/ConvertMQ0Bed.d \
./src/convert/Convert_Pindel.d \
./src/convert/Convert_VCF_to_BED.d \
./src/convert/Process_Lumpy.d 


# Each subdirectory must supply rules for building sources it contributes
src/convert/%.o: ../src/convert/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


