# This file is part of QCxMS-2.
set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs

"${dir}/xtb_mctc_accuracy.f90"
"${dir}/xtb_mctc_constants.f90"
"${dir}/xtb_mctc_convert.f90"
"${dir}/xtb_mctc_symbols.f90"
"${dir}/xtb_mctc_resize.f90"
"${dir}/boxmuller.f90"
"${dir}/iee.f90"
"${dir}/mcsimu.f90"
"${dir}/tsmod.f90"
"${dir}/iomod.f90"
"${dir}/qcxms2_data.f90"
"${dir}/argparser.f90"
"${dir}/fragmentation.f90"
"${dir}/qmmod.f90"
"${dir}/charges.f90"
"${dir}/utility.f90"
"${dir}/plot.f90"
"${dir}/readl.f90"
"${dir}/timing.f90"
"${dir}/reaction.f90"
"${dir}/isotopes.f90"
"${dir}/header.f90"
"${dir}/structools.f90"
"${dir}/rmsd.f90"
)
list(APPEND prog
  "${dir}/main.f90"
)
set(srcs ${srcs} PARENT_SCOPE)
set(prog ${prog} PARENT_SCOPE)

