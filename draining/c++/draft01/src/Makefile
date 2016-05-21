include ../Makefile.in

#==============================================================================
# all object files that must be linked into final executable
OBJS = 
OBJS += run.o 

#==============================================================================
# rule for compiling everything
all : run

#==============================================================================
# rule for building executable from object files
run : $(OBJS) $(LIBRARIES)
	$(LOADER) $(LOADOPTS) $(OBJS) \
	$(LIBRARIES) -o $@ -lstdc++
	./$@

#==============================================================================
# rule for compiling individual source files into object files
.cpp.o :
	$(COMPILE) -o $@ $<

#==============================================================================
# rule to clean up all output files
clean :
	rm -f *.o run
