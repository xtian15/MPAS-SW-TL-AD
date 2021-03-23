FC=f2py
LIB=module_sw_mpas.so

FLAG="-fdefault-real-8" # gnu95
#FLAG="-real-size 64" # intel

${LIB}:		mpas_sw_module.f90 mpas_sw_testcases.f90 mpas_sw_dyncore.f90 mpas_sw_operators.f90 mpas_sw_dyncore_tlm.f90 mpas_sw_dyncore_adj.f90 grid2mpas.so
		${FC} --fcompiler=gnu95 --f90flags=${FLAG} -m module_sw_mpas -c mpas_sw_module.f90 mpas_sw_dyncore.f90 mpas_sw_operators.f90 mpas_sw_dyncore_tlm.f90 mpas_sw_dyncore_adj.f90

grid2mpas.so:	grid2mpas.f90
		${FC} -m grid2mpas -c grid2mpas.f90
clean:
		rm -r ${LIB} ${LIB}.dSYM grid2mpas.so
pyclean:
		rm *.pyc
