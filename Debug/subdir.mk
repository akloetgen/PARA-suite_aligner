################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../QSufSort.c \
../bamlite.c \
../bntseq.c \
../bwa.c \
../bwamem.c \
../bwamem_pair.c \
../bwape.c \
../bwase.c \
../bwaseqio.c \
../bwt.c \
../bwt_gen.c \
../bwt_lite.c \
../bwtaln.c \
../bwtgap.c \
../bwtindex.c \
../bwtsw2_aux.c \
../bwtsw2_chain.c \
../bwtsw2_core.c \
../bwtsw2_main.c \
../bwtsw2_pair.c \
../example.c \
../fastmap.c \
../is.c \
../kopen.c \
../kstring.c \
../ksw.c \
../kthread.c \
../main.c \
../malloc_wrap.c \
../pemerge.c \
../utils.c 

OBJS += \
./QSufSort.o \
./bamlite.o \
./bntseq.o \
./bwa.o \
./bwamem.o \
./bwamem_pair.o \
./bwape.o \
./bwase.o \
./bwaseqio.o \
./bwt.o \
./bwt_gen.o \
./bwt_lite.o \
./bwtaln.o \
./bwtgap.o \
./bwtindex.o \
./bwtsw2_aux.o \
./bwtsw2_chain.o \
./bwtsw2_core.o \
./bwtsw2_main.o \
./bwtsw2_pair.o \
./example.o \
./fastmap.o \
./is.o \
./kopen.o \
./kstring.o \
./ksw.o \
./kthread.o \
./main.o \
./malloc_wrap.o \
./pemerge.o \
./utils.o 

C_DEPS += \
./QSufSort.d \
./bamlite.d \
./bntseq.d \
./bwa.d \
./bwamem.d \
./bwamem_pair.d \
./bwape.d \
./bwase.d \
./bwaseqio.d \
./bwt.d \
./bwt_gen.d \
./bwt_lite.d \
./bwtaln.d \
./bwtgap.d \
./bwtindex.d \
./bwtsw2_aux.d \
./bwtsw2_chain.d \
./bwtsw2_core.d \
./bwtsw2_main.d \
./bwtsw2_pair.d \
./example.d \
./fastmap.d \
./is.d \
./kopen.d \
./kstring.d \
./ksw.d \
./kthread.d \
./main.d \
./malloc_wrap.d \
./pemerge.d \
./utils.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


